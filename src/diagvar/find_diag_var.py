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
PRIMER_DELIM =  ('(', ')')
CRRNA_DELIM =  ('{', '}')
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
    # Get genotypes for each group
    group_dt = _group_genotypes(variant, groups)
    diag_gt = _diagnostic_gts(variant, groups)
        
    # Ignore sites that are in low mapping quality regions
    if variant.info['MQ'] < min_map_qual:
        return "Low mapping"
        
    
    # Ignore sites with not enough reads/samples
    group_map = {sample: group for group, samples in groups.items() \
                 for sample in samples}
    samples_count = Counter([group_map[x.name] for x in variant.samples.values() \
                            if x.name in group_map and x['DP'] >= min_reads and x['GQ'] >= min_geno_qual])
    not_enough = {g: g in samples_count and samples_count[g] < min_samples for g in groups}
    if any(not_enough.values()):
        return "Missing data"
        
    # Ignore spanning deletions
    diag_gt = {k: None if v == "*" else v for k, v in diag_gt.items()}

    # Ignore indels
    diag_gt = {k: v if v is not None and len(v) == len(variant.ref) else None for k, v in diag_gt.items()}

    # Ignore sites without enough canidate variants
    if sum([x != None for x in diag_gt.values()]) < min_groups:
        return "Few groups"
    
    # Return which variants are diagnostic for each group
    return diag_gt


class GroupedVariant:
    
    def __init__(self, variant, groups):
        self.variant = variant
        self.groups = groups
        self._counts = _group_genotypes(variant, groups)
        
    @staticmethod
    def _count_genotypes(variant, subset=None, hetero=True):
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
        if subset is None:
            subset = variant.samples.keys()
        
        counts = {}
        for sample_id, data in variant.samples.items():
            if subset is not None and sample_id not in subset:
                continue
            if data['DP'] == 0: # https://gatk.broadinstitute.org/hc/en-us/articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot
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
        return counts


    def _group_genotypes(variant, groups, hetero=True):
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
        return {group_name: _count_genotypes(variant, samples, hetero) \
                for group_name, samples in groups.items()}

    def _diagnostic_gts(variant, groups):
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
        # Get genotypes for each group
        geno_counts = _group_genotypes(variant, groups, hetero=False)
        # Get alleles for each groups
        alleles = {}
        for g in groups.keys():
            alleles[g] = set(geno_counts[g].keys())
            if UNKNOWN_CHAR in alleles[g]:
                alleles[g].remove(UNKNOWN_CHAR)
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

    
    def counts(self, hetero = True):
        output = copy.deepcopy(self._counts)
        if not hetero:
            for group, counts in self._counts.items():
                for allele, count in counts.items():
                    if HETERO_DELIM in allele:
                        hetero_parts = allele.split(HETERO_DELIM)
                        # Add counts for each haplotype
                        for part in hetero_parts:
                            if part in output[group]:
                                output[group][part] += count
                            else:
                                output[group][part] = count
                        # Remove the heterozygous count
                        del output[group][allele]
        return output


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
    object?
        information on regions that contain variants diagnostic for a
        subset of samples. TBD.
    """
    pass
