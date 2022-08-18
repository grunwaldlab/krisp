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


def _read_group_data(metadata_path, sample_col="sample_id", group_col="group"):
  """
  Reads metadata file and returns dictionary with group IDs as keys and lists of sample names as items
  """
  metadata = pandas.read_csv(metadata_path, sep='\t')
  output = {}
  for index, row in metadata.iterrows():
      group = row[group_col]
      sample = row[sample_col]
      if group in output:
          output[group].append(sample)
      else:
          output[group] = [sample]
  return output



def _check_variant(
        variant,
        groups,
        min_groups = 1,
        min_samples = 5,
        min_reads = 10,
        min_geno_qual = 40,
        min_map_qual = 50):
    """Check how a variant relates to a given subset of samples
  
    Classifies a variant into one of the following categories:
      * Diagnostic: conserved in the subset and different in the other samples
      * Conserved: conserved in the subset
      * Unconserved: not conserved in the subset
      * Missing data: Too few samples or reads to make a determination (`min_samples`, `min_reads`, `min_geno_qual`)
      * Low mapping: Too low mapping quality score (`min_map_qual`)
      * Few groups: The variant distinguishes too few of the groups (`min_groups`)
    
    Parameters
    ----------
    variant : pysam.libcbcf.VariantRecord
        The variant to process.
    targets : list of str
        The sample names representing the target group.
    nontargets : list of str
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
                 min_samples=5,
                 min_reads=10,
                 min_geno_qual=40):
        self.variant = variant
        self.groups = groups
        # self.groups = {g: [x for x in ids if x in variant.samples.keys()] for g, ids in groups.items()}
        self.min_samples = min_samples
        self.min_reads = min_reads
        self.min_geno_qual = min_geno_qual

        # Store counts of samples for each group
        #   Note: this can be different from the sum of counts in allele_counts below due to heterozygous positions
        self.sample_counts = self._sample_counts(variant, groups, min_reads=min_reads, min_geno_qual=min_geno_qual)

        # Store counts of each allele for each group
        self.allele_counts = self._allele_counts(variant, groups, min_reads=min_reads, min_geno_qual=min_geno_qual)

        # Store which alleles are conserved for each group
        self.conserved = self._conserved(variant, groups, min_samples=min_samples, min_reads=min_reads, min_geno_qual=min_geno_qual)

        # Store which alleles are diagnostic for each group
        self.diagnostic = self._diagnostic(variant, groups, min_samples=min_samples, min_reads=min_reads, min_geno_qual=min_geno_qual)

    @classmethod
    def _count_genotypes(cls, variant, subset=None, hetero=True, unknown=True, min_reads=0,
                         min_geno_qual=0):
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

        Returns
        -------
        dict of int
            The number of each genotype
        """
        # If no subset is specified, use all samples
        if subset is None:
            subset = variant.samples.keys()

        # Filter samples based on data quality
        read_counts = GroupedVariant._sample_read_counts(variant)
        geno_quals = GroupedVariant._sample_geno_qual(variant)
        subset = [s for s in subset \
                  if s in variant.samples.keys() \
                  if read_counts[s] >= min_reads \
                  and geno_quals[s] >= min_geno_qual]
        # Count alleles
        counts = {}
        for sample_id, data in variant.samples.items():
            if subset is not None and sample_id not in subset:
                continue
            if data[
                'DP'] == 0:  # https://gatk.broadinstitute.org/hc/en-us/articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot
                alleles = UNKNOWN_CHAR
            else:
                alleles = sorted(list(set(data.alleles)))
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
    def _allele_counts(cls, variant, groups, hetero=True, unknown=True, min_reads=10,
                       min_geno_qual=40):
        """For a given variant return the counts of each genotype for each group.

        Parameters
        ----------
        variant : pysam.libcbcf.VariantRecord
            The variant to process
        groups : dict of list of str
            The sample IDs within each group, named by group name.
        hetero : bool
            If `False`, heterozygous variants are counted once for each
            haplotype allele instead of counting once as a pair (e.g. "A/T")

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
                                                            min_geno_qual=min_geno_qual)
        return output

    @classmethod
    def _conserved(cls,
                   variant,
                   groups,
                   unknown=False,
                   min_samples=5,
                   min_reads=10,
                   min_geno_qual=40):
        """Check if an allele is conserved for each group

        Return
        ------
        dict of str/None
            Groups as keys, allele/None as value
        """
        # Get allele counts for each group
        geno_counts = cls._allele_counts(variant,
                                         groups,
                                         hetero=False,
                                         unknown=unknown,
                                         min_reads=min_reads,
                                         min_geno_qual=min_geno_qual)
        # Check if there is only one allele and there is enough reads
        samp_counts = cls._sample_counts(variant,
                                         groups,
                                         min_reads=min_reads,
                                        min_geno_qual=min_geno_qual)
        output = {}
        for group, counts in geno_counts.items():
            if len(counts) == 1 and samp_counts[group] >= min_samples:
                output[group] = list(counts.keys())[0]
            else:
                output[group] = None
        return output

    @classmethod
    def _diagnostic(cls, variant, groups, conserved = True,
                    min_samples = 5,
                    min_reads = 10,
                    min_geno_qual = 40):
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
        # If there are not enough samples for any group, no diagnostic variants can be found
        samp_counts = cls._sample_counts(variant, groups, min_reads=min_reads, min_geno_qual=min_geno_qual)
        if any([n < min_samples for n in samp_counts.values()]):
            return {group: None for group in groups.keys()}
        # Get genotypes for each group
        geno_counts = cls._allele_counts(variant, groups, hetero=False, unknown=False, min_reads=min_reads, min_geno_qual=min_geno_qual)
        # Get alleles for each groups
        alleles = {}
        for g in groups.keys():
            alleles[g] = set(geno_counts[g].keys())
            # if UNKNOWN_CHAR in alleles[g]:
            #     alleles[g].remove(UNKNOWN_CHAR)
        # Remove alleles for each group that appear in other groups
        diag = copy.deepcopy(alleles)
        for group in groups.keys():
            for other_group in groups.keys():
                if other_group != group:
                    diag[group] -= alleles[other_group]
        # Only return values for conserved diagnostic alleles
        for group in groups.keys():
            if len(alleles[group]) > 1 or len(diag[group]) == 0:
                diag[group] = None
            else:
                diag[group] = list(diag[group])[0]
        return diag
    
    @staticmethod
    def _sample_read_counts(variant):
        """Count the number of reads for each sample.
        
        Return
        ------
        dict of number
            The number of reads in each sample
        """
        return {s.name: s['DP'] for s in variant.samples.values()}
    
    @staticmethod
    def _sample_geno_qual(variant):
        """The genotype quality score for each sample.
        
        Return
        ------
        dict of number
            The genotype quality score in each sample
        """
        return {s.name: s['GQ'] for s in variant.samples.values()}
    
    @staticmethod
    def _subset_sample_counts(variant, subset, min_reads = 10, min_geno_qual = 40):
        """Number of samples passing filters in a given subset"""
        read_counts = GroupedVariant._sample_read_counts(variant)
        geno_quals = GroupedVariant._sample_geno_qual(variant)
        return(sum([read_counts[s] >= min_reads \
                    and geno_quals[s] >= min_geno_qual \
                    for s in subset if s in read_counts]))

    @classmethod
    def _sample_counts(cls, variant, groups, min_reads = 10, min_geno_qual = 40):
        """Number of samples in each group"""
        output = {}
        for group, samples in groups.items():
            output[group] = cls._subset_sample_counts(variant, samples,
                                                       min_reads=min_reads, min_geno_qual=min_geno_qual)
        return output


    def all_conserved(self,
                       unknown = False,
                       min_samples = 5,
                       min_reads = 10,
                       min_geno_qual = 40):
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
        nontarget = None,
        min_groups = 1,
        min_samples = 5,
        min_reads = 10,
        min_geno_qual = 40,
        min_map_qual = 50):
    """Find regions with diagnostic variants
    
    Return information on regions that contain variants diagnostic for a subset
    of samples, optionally surrounded by conserved regions where primers can be
    designed.
    
    Parameters
    ----------
    variants : an iterable returning variants or str
        A series of consecutive variants in which to identifiy diagnostic 
        clusters of variants.
    groups : list/dict/tuple of list/tuple of str
        The sample IDs for each group to find diagnostic clusters for.
    nontarget : list of str, optional
        The sample IDs that should be distict from those in `groups`, but
        are otherwise not of interest. `min_samples` does not apply to these.
    min_groups : int, optional
        The minimum number of `groups` that must be uniquly distinguished by
        each variants. Increasing this will significantly reduce results.
    min_samples : int, optional
        The minimum number of samples that must represent each group in
        `groups` for a given variants. Samples must pass the `min_reads` 
        and `min_geno_qual` filters to count twoards this minimum.
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

