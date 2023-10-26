import sys
import argparse
import pandas
import gzip
import copy
import primer3
import pysam
import re
from collections import deque
from Bio import SeqIO
from Bio import Seq
from Bio.Data import IUPACData
from contextlib import contextmanager
from collections import Counter

# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM = ('(', ')')
CRRNA_DELIM = ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'


def _parse_group_data(metadata_path, groups=None, sample_col="sample_id", group_col="group", min_samples=None):
    """Reads metadata file and returns dictionary with group IDs as keys and
    lists of sample names as items
    """
    # Read samples in each group
    metadata = pandas.read_csv(metadata_path, sep=',')
    output = {}
    for index, row in metadata.iterrows():
        group = row[group_col]
        sample = row[sample_col]
        if group in output:
            output[group].append(sample)
        else:
            output[group] = [sample]
    # Check if user-defined groups are present
    if groups is not None:
        missing_groups = [g for g in groups if g not in output.keys()]
        if len(missing_groups) > 0:
            raise ValueError((
                f'One or more user-defined groups are not present in the metadata file:\n'
                f'    {metadata_path}\n'
                f'The following user-defined groups are not present:\n'                
                f'    {", ".join(missing_groups)}\n'
                f'The following groups are present in the metadata file:\n'
                f'    {", ".join(output.keys())}'
            ))
    # Check that there are enough samples in all the groups
    if min_samples is not None:
        groups_with_too_few_samples = {k: len(v) for k, v in output.items() if k in groups and len(v) < min_samples}
        if len(groups_with_too_few_samples) > 0:
            raise ValueError((
                f'One or more user-defined groups have fewer samples than `--min_samples`:\n'
                f'    {", ".join([g + " (" + str(c) + ")" for g, c in groups_with_too_few_samples.items()])}'
            ))
    # Subset to just groups of interest
    if groups is not None:
        output = {g: v for g, v in output.items() if g in groups}
    return output


def _check_variant(
        variant,
        groups,
        targets,
        min_groups=1,
        min_samples=5,
        min_reads=10,
        min_geno_qual=40,
        min_map_qual=50):
    """Check how a variant relates to a given subset of samples

    Classifies a variant into one of the following categories:
      * Diagnostic: conserved in the subset and different in the other samples
      * Conserved: conserved in the subset
      * Unconserved: not conserved in the subset
      * Missing data: Too few samples or reads to make a determination
      (`min_samples`, `min_reads`, `min_geno_qual`)
      * Low mapping: Too low mapping quality score (`min_map_qual`)
      * Few groups: The variant distinguishes too few of the groups
      (`min_groups`)

    Parameters
    ----------
    variant : pysam.libcbcf.VariantRecord
        The variant to process.
    targets : list of str
        The sample names representing the target group.
    nontargets: list of str
        The sample names representing the group to be distinguished
        from the target.
    """
    pass
    # # Get genotypes for each group
    # group_dt = _group_genotypes(variant, groups)
    # diag_gt = _diagnostic_gts(variant, groups)
    #
    # # Ignore sites that are in low mapping quality regions
    # if variant.info['MQ'] < min_map_qual:
    #     return "Low mapping"
    #
    #
    # # Ignore sites with not enough reads/samples
    # group_map = {sample: group for group, samples in groups.items() \
    #              for sample in samples}
    # samples_count = Counter([group_map[x.name] for x in variant.samples.values() \
    #                         if x.name in group_map and x['DP'] >= min_reads and x['GQ'] >= min_geno_qual])
    # not_enough = {g: g in samples_count and samples_count[g] < min_samples for g in groups}
    # if any(not_enough.values()):
    #     return "Missing data"
    #
    # # Ignore spanning deletions
    # diag_gt = {k: None if v == "*" else v for k, v in diag_gt.items()}
    #
    # # Ignore indels
    # diag_gt = {k: v if v is not None and len(v) == len(variant.ref) else None for k, v in diag_gt.items()}
    #
    # # Ignore sites without enough canidate variants
    # if sum([x != None for x in diag_gt.values()]) < min_groups:
    #     return "Few groups"
    #
    # # Return which variants are diagnostic for each group
    # return diag_gt


class GroupedVariant:
    
    def __init__(self, variant, groups,
                 check_groups=False,
                 min_samp_prop=0.9,
                 min_samples=5,
                 min_reads=10,
                 min_geno_qual=40,
                 min_freq=0.1,
                 min_map_qual=30,
                 min_var_qual=10,
                 force=False):
        self.variant = variant
        if check_groups:  # optional check needed for reading from stdin
            metadata_samples = {val for values in groups.values() for val in values}
            vcf_samples = set(variant.samples.keys())
            missing_in_vcf = metadata_samples - vcf_samples
            missing_in_meta = vcf_samples - metadata_samples
            if len(missing_in_vcf) > 0 and not force:
                raise ValueError(f'The following samples specified in the metadata cannot be found in the VCF input:\n'
                                 f'    {", ".join(missing_in_vcf)}\n')
            self.groups = {g: [x for x in ids if x in variant.samples.keys()] for g, ids in groups.items()}
        else:
            self.groups = groups
        self.min_samples = min_samples
        self.min_reads = min_reads
        self.min_geno_qual = min_geno_qual
        self.min_freq = min_freq

        # Store basic info

        # Store counts of samples for each group
        #   Note: this can be different from the sum of counts in allele_counts
        #   below due to heterozygous positions
        count_data = self._sample_counts(variant, self.groups,
                                                 min_reads=min_reads,
                                                 min_geno_qual=min_geno_qual)
        self.sample_counts = count_data['counts']
        self.missing_samp_ids = count_data['missing']

        # Store counts of each allele for each group
        self.allele_counts = self._allele_counts(variant, self.groups,
                                                 hetero=False,
                                                 min_reads=min_reads,
                                                 min_geno_qual=min_geno_qual,
                                                 min_freq=min_freq)

        # Store which alleles are conserved for each group
        self.conserved = self._conserved(min_samp_prop=min_samp_prop,
                                         min_samples=min_samples,
                                         min_map_qual=min_map_qual,
                                         min_var_qual=min_var_qual)

        # Store which alleles are diagnostic for each group
        self.diagnostic = self._diagnostic(min_samp_prop=min_samp_prop,
                                           min_samples=min_samples,
                                           min_map_qual=min_map_qual,
                                           min_var_qual=min_var_qual)

    @classmethod
    def from_vcf(cls, variants, groups, **kwargs):
        """Iterate over a variants, creating GroupedVariants objects

        Groups are checked once per vcf file to enhance performance
        """
        groups_checked = False
        for var in variants:
            if groups_checked:
                out = cls(var, groups, check_groups=False, **kwargs)
            else:
                out = cls(var, groups, check_groups=True, **kwargs)
                groups = out.groups
                groups_checked = True
            yield out

    @classmethod
    def _count_genotypes(cls, variant, subset=None, hetero=True, unknown=True,
                         min_reads=0, min_geno_qual=0, min_freq=0.1):
        """For a variant return the counts of genotypes for a subset of samples.

        Parameters
        ----------
        variant : pysam.libcbcf.VariantRecord
            The variant to process
        subset : list of str, optional
            The sample names to count genotypes for. If `None`, use all samples.
        hetero : bool
            If `False`, heterozygous variants are counted once for each
            haplotype allele instead of counting once as a pair (e.g. "A/T")
        min_reads : int, optional
            The minimum number of reads a sample have at the location of a
            given variant.
        min_geno_qual : int, optional
            The minimum genotype quality score (phred scale). This corresponds
            to the per-sample output in the VCF encoded as "GQ".
        min_freq : float, optional
            The minimum proportion of reads an allele must have to be considered real.
            This is meant to counter sequencing errors.
            If not set, the genotypes called in the VCF are used without considering read depth.

        Returns
        -------
        dict of int
            The number of each genotype
        """
        # If no subset is specified, use all samples
        if subset is None:
            subset = variant.samples.keys()

        # Filter samples based on data quality
        subset = [s for s in subset
                  if variant.samples[s]['DP'] is not None
                  and variant.samples[s]['DP'] >= min_reads
                  and variant.samples[s]['GQ'] is not None
                  and variant.samples[s]['GQ'] >= min_geno_qual]

        # Count alleles
        counts = {}
        for sample_id, data in variant.samples.items():
            if subset is not None and sample_id not in subset:
                continue
            if data['DP'] == 0:  # https://gatk.broadinstitute.org/hc/en-us/
                # articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot
                alleles = UNKNOWN_CHAR
            else:
                if min_freq is None:
                    alleles = sorted(list(set(data.alleles)))
                else:
                    min_depth = sum(data['AD']) * min_freq
                    alleles = sorted(list(set([variant.alleles[i] for i, d in enumerate(data['AD']) if d > 0 and d >= min_depth])))
                alleles = [UNKNOWN_CHAR if a is None else a for a in alleles]
                if hetero:
                    alleles = ["/".join(alleles)]
            for allele in alleles:
                if allele in counts:
                    counts[allele] += 1
                else:
                    counts[allele] = 1
        # Remove ambiguous alleles
        if not unknown:
            counts = {k: v for k, v in counts.items() if k != "?"}
        return counts

    @classmethod
    def _allele_counts(cls, variant, groups, hetero=True, unknown=True,
                       min_reads=10, min_geno_qual=40, min_freq=0.1):
        """For a given variant return the counts of each genotype for each
        group.

        Parameters
        ----------
        variant : pysam.libcbcf.VariantRecord
            The variant to process
        groups : dict of list of str
            The sample IDs within each group, named by group name.
        hetero : bool
            If `False`, heterozygous variants are counted once for each
            haplotype allele instead of counting once as a pair (e.g. "A/T")
        min_reads : int, optional
            The minimum number of reads a sample have at the location of a
            given variant.
        min_geno_qual : int, optional
            The minimum genotype quality score (phred scale). This corresponds
            to the per-sample output in the VCF encoded as "GQ".
        min_freq : float, optional
            The minimum proportion of reads an allele must have to be considered real.
            This is meant to counter sequencing errors.
            If not set, the genotypes called in the VCF are used without considering read depth.

        Returns
        -------
        dict of dict of int
            For each group: the number of counts for each allele.
        """
        output = {}
        for group, samples in groups.items():
            output[group] = GroupedVariant._count_genotypes(variant,
                                                            samples,
                                                            hetero=hetero,
                                                            unknown=unknown,
                                                            min_reads=min_reads,
                                                            min_geno_qual=min_geno_qual,
                                                            min_freq=min_freq)
        return output

    def _conserved(self, min_samp_prop=0.9, min_samples=5, min_map_qual=30, min_var_qual=10):
        """Check if an allele is conserved for each group

        Return
        ------
        dict of str/None
            Groups as keys, allele/None as value
        """
        # Cant determine if RMS mapping quality is too low
        if self.variant.info['MQ'] < min_map_qual:
            return {group: None for group in self.groups.keys()}

        # Cant determine if variant quality (QUAL column) is too low
        if self.variant.qual < min_var_qual:
            return {group: None for group in self.groups.keys()}

        # Cant determine if conserved if too few samples
        output = {}
        for group, counts in self.allele_counts.items():
            samp_prop = self.sample_counts[group] / len(self.groups[group])
            if len(counts) == 1 and self.sample_counts[group] >= min_samples and samp_prop >= min_samp_prop:
                output[group] = list(counts.keys())[0]
            else:
                output[group] = None
        return output

    def _diagnostic(self, min_samp_prop=0.9, min_samples=5, min_map_qual=30, min_var_qual=10):
        """
        Get conserved variants only present in each group.
        
        Parameters
        ----------
        variant : pysam.libcbcf.VariantRecord
            The variant to process
        groups : dict of list of str
            The sample IDs within each group, named by group name.
        
        Returns
        -------
        dict of str
            For each group, the alleles that are conserved and unique to
            that group.
        """
        # Cant determine if RMS mapping quality is too low
        if self.variant.info['MQ'] < min_map_qual:
            return {group: None for group in self.groups.keys()}

        # Cant determine if variant quality (QUAL column) is too low
        if self.variant.qual < min_var_qual:
            return {group: None for group in self.groups.keys()}

        # If there are not enough samples for any group, no diagnostic variants found
        if any([n < min_samples or n / len(self.groups[g]) < min_samp_prop for g, n in self.sample_counts.items()]):
            return {group: None for group in self.groups.keys()}

        # Get alleles for each groups
        alleles = {}
        for g in self.groups.keys():
            alleles[g] = set(self.allele_counts[g].keys())
            # if UNKNOWN_CHAR in alleles[g]:
            #     alleles[g].remove(UNKNOWN_CHAR)
        # Remove alleles for each group that appear in other groups
        diag = copy.deepcopy(alleles)
        for group in self.groups.keys():
            for other_group in self.groups.keys():
                if other_group != group:
                    diag[group] -= alleles[other_group]
        # Only return values for conserved diagnostic alleles
        for group in self.groups.keys():
            if len(alleles[group]) > 1 or len(diag[group]) == 0:
                diag[group] = None
            else:
                diag[group] = list(diag[group])[0]
        return diag

    @staticmethod
    def _subset_sample_counts(variant, subset, min_reads=10, min_geno_qual=40):
        """Number of samples passing filters in a given subset"""
        is_good = {s: variant.samples[s]['DP'] is not None
                    and variant.samples[s]['DP'] >= min_reads
                    and variant.samples[s]['GQ'] is not None
                    and variant.samples[s]['GQ'] >= min_geno_qual
                    for s in subset}
        missing_samp_ids = {k for k, v in is_good.items() if not v}
        return {"counts": sum(is_good.values()), "missing": missing_samp_ids}

    @classmethod
    def _sample_counts(cls, variant, groups, min_reads=10, min_geno_qual=40):
        """Number of samples in each group"""
        counts = {}
        missing_samp_ids = {}
        for group, samples in groups.items():
            output = cls._subset_sample_counts(variant, samples,
                                                      min_reads=min_reads,
                                                      min_geno_qual= min_geno_qual)
            counts[group] = output['counts']
            missing_samp_ids[group] = output['missing']
        return {'counts': counts, 'missing': missing_samp_ids}

    def allele_lens(self, group):
        """Number of nucleotides of each allele

        Returns
        -------
        dict of int
            The number of nucleotides for each allele, named by allele
        """
        alleles = self.allele_counts[group].keys()
        out = {}
        for allele in alleles:
            if "/" in allele:
                out[allele] = max([len(x) for x in allele.split("/")])
                continue
            if allele == "*":
                out[allele] = 0
                continue
            out[allele] = len(allele)
        return out

    def max_allele_len(self, group):
        """Maximum length of an allele for a given group"""
        if len(self.allele_counts[group]) == 0:
            return len(self.variant.ref)
        else:
            return max(self.allele_lens(group).values())

    def all_conserved(self,
                      unknown=False,
                      min_samples=5,
                      min_reads=10,
                      min_geno_qual=40):
        """Check if an allele is conserved in all groups
        
        Return
        ------
        str/None
            If all groups have the same variant, return the allele,
            otherwise return None
        """
        pass
        # is_consv = self.conserved(unknown = unknown,
        #                           min_samples = min_samples,
        #                           min_reads = min_reads,
        #                           min_geno_qual = min_geno_qual)
        # alleles = set(is_consv.values())
        # if len(alleles) == 1:
        #     return list(alleles)[0]
        # else:
        #     return None


def find_diag_var(
        variants,
        groups,
        nontarget=None,
        min_groups=1,
        min_samples=5,
        min_reads=10,
        min_geno_qual=40,
        min_map_qual=50):
    """Find regions with diagnostic variants
    
    Return information on regions that contain variants diagnostic for a subset
    of samples, optionally surrounded by conserved regions where primers can be
    designed.
    
    Parameters
    ----------
    variants : an iterable returning variants or str
        A series of consecutive variants in which to identify diagnostic
        clusters of variants.
    groups : list/dict/tuple of list/tuple of str
        The sample IDs for each group to find diagnostic clusters for.
    nontarget : list of str, optional
        The sample IDs that should be distinct from those in `groups`, but
        are otherwise not of interest. `min_samples` does not apply to these.
    min_groups : int, optional
        The minimum number of `groups` that must be uniquely distinguished by
        each variant. Increasing this will significantly reduce results.
    min_samples : int, optional
        The minimum number of samples that must represent each group in
        `groups` for a given variants. Samples must pass the `min_reads` 
        and `min_geno_qual` filters to count towards this minimum.
    min_reads : int, optional
        The minimum number of reads a sample have at the location of a
        given variant.
    min_geno_qual : int, optional
        The minimum genotype quality score (phred scale). This corresponds
        to the per-sample GATK output in the VCF encoded as "GQ".
    min_map_qual : int, optional
        The minimum root square mean mapping quality of reads supporting a
        variant (phred scale). This corresponds to the per-variant GATK 
        output in the VCF encoded as "MQ".
        
    Yields
    ------
    GroupedVariant objects for variants that pass the filters 
    """
    pass
