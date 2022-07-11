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


def classify_variant(record, subset, others):
    """Check how a variant relates to a given subset of samples
  
    Classifies a variant into one of the following categories:
      * Diagnostic: conserved in the subset and different in the other samples
      * Conserved: conserved in the subset
      * Unconserved: not conserved in the subset
      * Missing data: Too few samples or reads to make a determination
      * Low mapping: Too low mapping quality score
    """


def find_diag_var(
        variants,
        groups,
        nontarget = None,
        min_groups = 1,
        min_samples = 5,
        min_reads = 10,
        min_geno_qual = 40,
        min_map_qual = 50,
        min_freq = 0.95):
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
    min_freq : float, optional
        The proportion of reads a variant must have to be considered real.
        This is meant to counter sequencing errors.
        
    Yields
    ------
    object?
        information on regions that contain variants diagnostic for a
        subset of samples. TBD.
    """
