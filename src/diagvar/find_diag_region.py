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

from .find_diag_var import find_diag_var

# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM =  ('(', ')')
CRRNA_DELIM =  ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'

class var_sliding_window:
    """Stores all the queues so that each group can have its own spacer queue"""
    def __init__(self, sample_subset, spacer_len, flank_len):
        self.sample_subset = sample_subset
        self.spacer_len = spacer_len
        self.flank_len = flank_len
        self.upstream = deque(maxlen=flank_len)
        self.downstream = deque(maxlen=flank_len)
        self.spacer = deque()
    
    def add_variant(self, variant):
        if len(self.upstream) < self.flank_len:
            # Pre-load upstream queue before processing any records
            self.upstream.append(variant)
        else:
            # Move one variant to spacer and add next to upstream queue
            self.spacer.append(self.upstream[0])
            self.upstream.append(variant)
            # Move from spacer to downstream queue once spacer is too large
            while self._var_span_len(self.spacer) > self.spacer_len:
                self.downstream.appendleft(self.spacer.pop(0))
    
    def _var_span_len(self, variants):
        # TODO: this will not handle indels correctly, needs to be subset-specific.
        positions = [var.pos for var in variants]
        return max(positions) - min(positions)


def _check_variant_cluster(variants, subset):
    """Performs a series of checks to see if the cluster is diagnostic"""
    
    # Is there conserved sequence for a spacer?
    
    # Are there enough diagnostic sites?
    
    # Can primers be designed in adjacent conserved regions?
    
    pass

def find_diag_region(
        variants,
        groups,
        nontarget = None,
        reference = None,
        primer3 = False,
        min_vars = 1,
        min_bases = 1,
        min_groups = 1,
        min_samples = 5,
        min_reads = 10,
        min_geno_qual = 40,
        min_map_qual = 50,
        min_freq = 0.95,
        spacer_len = 28,
        snp_offset = -3):
    """Find regions with diagnostic variants
    
    Return information on regions that contain variants diagnostic for a subset
    of samples, optionally surrounded by conserved regions where primers can be
    designed.
    
    Parameters
    ----------
    variants : an iterable returning variants or str
        A series of consecutive variants in which to identifiy diagnostic 
        clusters of variants.
    groups : dict of list/tuple of str
        The sample IDs for each group to find diagnostic clusters for,
        named by group.
    nontarget : list of str, optional
        The sample IDs that should be distict from those in `groups`, but
        are otherwise not of interest. `min_samples` does not apply to these.
    reference : Bio.SeqRecord.SeqRecord or str, optional
        A reference sequence to use or the path to a FASTA file with 
        a single sequence. Required if Primer3 is to be used.
    primer3 : bool, optional
        If `True`, run Primer3 on each potential cluster of diagnostic
        variants and only return results for which primers can be found
        in adjacent conserved sequence.
    min_vars : int, optional
        The minimum number of diagnostic variants. Note that a single
        indel variant can represent multiple diagnostic nucleotides.
        See `min_bases` to filter based on number of nucleotides.
    min_bases : int, optional
        The minimum number of diagnostic nucleotides. Note that a single
        indel variant can represent multiple diagnostic nucleotides.
        See `min_vars` to filter based on number of variants.
    min_groups : int, optional
        The minimum number of `groups` that must be uniquly distinguished by
        the same cluster of variants.
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
    spacer_len: int, optional
        The length of the spacer sequence that the crRNA will bind to and
        the primers will amplify.
    var_offset : int, optional
        How many bases from the 5` end of the spacer will the SNP occur.
        For example 0 is the 5` end of the spacer and -1 is the 3` end.
        
    Yields
    ------
    object?
        information on regions that contain variants diagnostic for a
        subset of samples. TBD.
    """
    flank = 100
    windows = {g: var_sliding_window(subset, spacer_len=spacer_len, flank_len=flank) for g, subset in groups.items()}
    for variant in find_diag_var():
        for group, subset in groups.items():
            windows[group].add_variant(variant)
            _check_variant_cluster(windows[group].spacer, subset)
            
        
     
