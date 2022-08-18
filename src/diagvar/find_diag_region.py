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

from .find_diag_var import find_diag_var, GroupedVariant

# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM =  ('(', ')')
CRRNA_DELIM =  ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'
iupac_key = {tuple((x for x in sorted(v))):k for k,v in IUPACData.ambiguous_dna_values.items()}
iupac_key[(UNKNOWN_CHAR, )] = 'N'

def collapse_to_iupac(seqs):
    """Combine sequences into a consensus using IUPAC ambiguity codes
    
    Parameters:
    -----------
    segs : list of str
        The sequences to combine.
    
    Returns
    -------
    str
        The combined sequence
    """
    seq_lens = [len(x) for x in seqs]
    max_len = max(seq_lens)
    if len(set(seq_lens)) != 1: #TODO: replace with alignment
        return '-' * max_len
    output = []
    for i in range(max_len):
        column = {s[i] for s in seqs}
        if "*" in column:
            output.append('N')
        else:
            output.append(iupac_key[tuple(sorted(column))])
    return "".join(output)
    

class GroupedRegion:
    
    def __init__(self, variants, groups):
        """
        Parameters
        ----------
        variants : list of 
        """
        #self.variants = [GroupedVariant(v, groups) for v in variants]
        self.variants = variants
        self.groups = groups
    
    @staticmethod
    def _get_reference(ref_path):
        if ref_path is None:
            return None
        else:
            if ref_path.endswith('.gz'):
                handle = gzip.open(ref_path, "rt")
            else:
                handle = open(ref_path)
                reference = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                names = [rec.id for rec in reference]
                reference = dict(zip(names, reference))
            return reference

    
    @staticmethod
    def region_length(variants, subset):
        """Get the sequence length spanned by variants in a subset
        
        Return
        ------
        int/None
            The number of nucelotides spanned by the variants. If any
            of the samples have variable length alleles, then a single
            length cannot be determined and None is returned.
        """
    
    @staticmethod
    def _count_genotypes(variants, subset, hetero=True, unknown=True):
        """For a sample subset, count each allele of each variant.
        
        Returns
        -------
        list of dict of int
            For each variant, the counts of each allele named by allele
        """
        output = []
        for variant in variants:
            output.append(GroupedVariant._count_genotypes(variant, subset=subset, hetero=hetero, unknown=unknown))
        return output
    
    @staticmethod
    def _subset_sequence(variants, subset,
                         reference = None,
                         counts = True,
                         hetero = True,
                         unknown = True,
                         min_reads = 0,
                         min_geno_qual = 0,
                         min_samples = 0,
                         consensus = False,
                         conserved = False,
                         wrap = '<>',
                         sep = ','):
        """For a sample subset, make a sequence representation.
        Parameters
        ----------
        reference : dict of str
            If supplied, non-variant reference positions are also
            returned. made by _get_reference
        counts
            If True, append sample counts to each allele
        hetero
            If False, heterozygous positions are represented by
            IUPAC ambiguity codes instead of being delimited with /
        consensus
            If True, combine all sample alleles into a single allele 
            using IUPAC codes. When a mixture of indel lengths is 
            present, TODO not sure, for now '-' * max indel length
        wrap
            If not None or '', a string composed of the first and last
            character used to border sets of multiple alleles. Has
            no effect if consensus is True.
        sep The character used to seperate multiple alleles.
        
        Returns
        -------
        list of str
            A representation of a sequence, with one element per
            variant or single base from the reference sequence.
        """
        # Create rendered text for variant information
        rendered = []
        for variant in variants:
            allele_counts = GroupedVariant._count_genotypes(variant, subset=subset, hetero=hetero, unknown=unknown,
                               min_reads=min_reads, min_geno_qual=min_geno_qual)
            sample_count = GroupedVariant._subset_sample_counts(variant, subset,
                                                      min_reads=min_reads,
                                                      min_geno_qual=min_geno_qual)
            if len(allele_counts) == 0:
                allele_counts = {'?':0}
            if conserved and (len(allele_counts) != 1 or sample_count < min_samples):
                text = None
            else:
                if consensus:
                    text = collapse_to_iupac(allele_counts.keys())
                else:
                    if counts:
                        text_parts = [a + str(c) for a, c in allele_counts.items()]
                    else:
                        text_parts = allele_counts.keys()
                    text = sep.join(text_parts)
                    if wrap is not None and wrap != '':
                        text = wrap[0] + text + wrap[1]
            rendered.append(text)
        # Insert into reference sequence
        if reference is None:
            return rendered
        else:
            # Get reference sequence spanning variants
            var_pos = [var.pos for var in variants]
            var_chrom = {var.chrom for var in variants}
            if len(var_chrom) > 1:
                raise ValueError('Variants cannot span multiple chromosomes')
            chrom = list(var_chrom)[0]
            ref_start = min(var_pos)
            ref_end = max(var_pos)
            ref_seq = list(reference[chrom].seq[ref_start:ref_end])
            # Replace reference with each variant
            for var, rend in zip(variants, rendered):
                replace_start = var.pos - ref_start
                replace_end = replace_start + len(var.ref)
                ref_seq = ref_seq[:replace_start] + [rend] + ref_seq[replace_end:]
            return ref_seq
    
    def sequence(self, reference = None, counts = False, ):
        """Infer the sequence for each group
        """
    
    def conserved(self, group = None, reference = None):
        """For each variant return alleles conserved in a given group
        
        Returns
        -------
        list of str/None
            Alleles for each variant when there is a conserved allele,
            otherwise None.
        """



class var_sliding_window:
    """Stores all the queues so that each group can have its own spacer queue"""
    def __init__(self, subset, spacer_len, flank_len, reference = None, min_reads = 0,
                         min_geno_qual = 0,
                         min_samples = 0):
        self.subset = subset
        self.spacer_len = spacer_len
        self.flank_len = flank_len
        self.upstream = deque()
        self.downstream = deque(maxlen=flank_len)
        self.spacer = deque()
        self.min_reads = min_reads
        self.min_geno_qual = min_geno_qual
        self.min_samples = min_samples
        self.reference = reference
    
    def increment(self):
        # Move one variant to spacer and add next to upstream queue
        self.spacer.append(self.upstream.popleft())
        # Move from spacer to downstream queue once spacer is too large
        while self._spacer_span() > self.spacer_len:
            self.downstream.appendleft(self.spacer.popleft())
    
    # def is_primed(self):
        # return len(self.upstream) >= self.flank_len
        
    # def add_variant(self, variant):
        # self.upstream.append(variant)
        # if self.is_primed():
            # self.increment()
        
    def _spacer_span(self):
        return len(''.join(self.spacer_seq()))


def window_generator(min_samples, min_reads, min_geno_qual, spacer_len, flank_len = 100):
    filter_args = {"min_samples": min_samples,
        "min_reads": min_reads,
        "min_geno_qual": min_geno_qual}
    # Initialize sliding windows
    windows = {}
    for group, subset in groups.items():
        windows[group] = var_sliding_window(subset, spacer_len=spacer_len,
                                            flank_len=flank, reference = ref, **filter_args)
    # Read through variants and put into windows
    for index, variant in enumerate(variants):
        for group in groups.keys():
            windows[group].upstream.append(variant)
            if index + 1 >= flank_len: # Once the upstream is full
                windows[group].increment()
                yield group, windows[group]
    # Clear the upstream window
    for index in range(flank_len - 1):
        for group in groups.keys():
            windows[group].increment()
            yield group, windows[group]


def _check_variant_cluster(variants, subset):
    """Performs a series of checks to see if the cluster is diagnostic"""
    
    
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
        snp_offset = 2):
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
        The minimum ngroupsumber of diagnostic nucleotides. Note that a single
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
    offset_left : int, optional
        The minimum number of bases on the 3' end of the spacer before
        the first diagnositic SNP.
    offset_right : int, optional
        The number of bases on the 5' end after the last diagnostic SNP.
        
    Yields
    ------
    object?
        information on regions that contain variants diagnostic for a
        subset of samples. TBD.
    """
    flank = 100 # TODO: base on max amplicon size
    window_width = spacer_len - offset_right - offset_left
    ref = GroupedRegion._get_reference('test_data/PR-102_v3.1.fasta')
    vcf_reader = GroupedVariant.from_vcf(variants, groups, min_samples=min_samples, min_reads=min_reads, min_geno_qual=min_geno_qual, min_map_qual=min_map_qual)
    windower = window_generator(vcf_reader, groups, window_width=window_width, flank_len=flank)
    for group, window in windower:
        region = GroupedRegion(window.spacer, groups)
        # Are all the variants in the spacer conserved?
        if any([x is None for x in region.conserved(group)]):
            continue
        # Is there conserved sequence upstream for the 5' end?
        nearby_vars = get_nearby_vars(window.upstream, offset_right) +\
                      get_nearby_vars(window.downstream, offset_left)
        if not all([v.is_group_conserved(group) for v in nearby_vars]):
            continue
        # Are there enough diagnostic variants?
        n_diag_var = sum([x is not None for x in region.diagnostic(group)])
        if n_diag_var < min_vars:
            continue
        # Can primers be designed in adjacent conserved regions?
        infer_adjacent_seq(window.upstream, max_len = max_amplicon_length)
        run_primer3(
    
    
                
            
            
        
     
