import logging
import os.path
import sys
import argparse
import gzip
import primer3
import pysam
import re
import itertools
import multiprocessing as mp
import queue

from collections import deque, defaultdict
from Bio import SeqIO
from Bio.Data import IUPACData
from contextlib import contextmanager
from logging.handlers import QueueHandler
from statistics import mean
from nltk.metrics import distance

from .find_diag_var import GroupedVariant, _parse_group_data
from .print_align import render_variant, Annotation

# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM = ('(', ')')
CRRNA_DELIM = ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'
iupac_key = {tuple((x for x in sorted(v))): k for k, v in
             IUPACData.ambiguous_dna_values.items()}
iupac_key[(UNKNOWN_CHAR,)] = 'N'
# primer3_col_names = [
#     'PRIMER_PAIR_0_PRODUCT_SIZE',
#     'PRIMER_PAIR_0_PENALTY',
#     'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 'PRIMER_INTERNAL_0_SEQUENCE',
#     'PRIMER_LEFT_0_PENALTY', 'PRIMER_RIGHT_0_PENALTY', 'PRIMER_INTERNAL_0_PENALTY',
#     'PRIMER_LEFT_0_TM', 'PRIMER_RIGHT_0_TM', 'PRIMER_INTERNAL_0_TM',
#     'PRIMER_LEFT_0_GC_PERCENT', 'PRIMER_RIGHT_0_GC_PERCENT', 'PRIMER_INTERNAL_0_GC_PERCENT',
#     'PRIMER_LEFT_0_SELF_ANY_TH', 'PRIMER_RIGHT_0_SELF_ANY_TH', 'PRIMER_INTERNAL_0_SELF_ANY_TH',
#     'PRIMER_LEFT_0_SELF_END_TH', 'PRIMER_RIGHT_0_SELF_END_TH', 'PRIMER_INTERNAL_0_SELF_END_TH',
#     'PRIMER_LEFT_0_HAIRPIN_TH', 'PRIMER_RIGHT_0_HAIRPIN_TH', 'PRIMER_INTERNAL_0_HAIRPIN_TH',
#     'PRIMER_LEFT_0_END_STABILITY', 'PRIMER_RIGHT_0_END_STABILITY',
#     'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH',
# ]
primer3_col_names = [
    'PRIMER_PAIR_0_PRODUCT_SIZE',
    'PRIMER_PAIR_0_PENALTY',
    'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE',
    'PRIMER_LEFT_0_PENALTY', 'PRIMER_RIGHT_0_PENALTY',
    'PRIMER_LEFT_0_TM', 'PRIMER_RIGHT_0_TM',
    'PRIMER_LEFT_0_GC_PERCENT', 'PRIMER_RIGHT_0_GC_PERCENT',
    'PRIMER_LEFT_0_SELF_ANY_TH', 'PRIMER_RIGHT_0_SELF_ANY_TH',
    'PRIMER_LEFT_0_SELF_END_TH', 'PRIMER_RIGHT_0_SELF_END_TH',
    'PRIMER_LEFT_0_HAIRPIN_TH', 'PRIMER_RIGHT_0_HAIRPIN_TH',
    'PRIMER_LEFT_0_END_STABILITY', 'PRIMER_RIGHT_0_END_STABILITY',
    'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH',
]
primer3_col_key = {n: n.replace("PRIMER_", "").replace("_0", "").lower() for n in primer3_col_names}

# Globals
logger = logging.getLogger(__name__)
failure_event = None

def configure_global_logger(args=None, mode="w"):
    global logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    stderr_handler = logging.StreamHandler()
    stderr_handler.setLevel(logging.WARNING)
    log_formatter = logging.Formatter('%(levelname)s: %(name)s: %(message)s')
    stderr_handler.setFormatter(log_formatter)
    logger.addHandler(stderr_handler)
    if args is not None:
        if args.log is None:
            if args.log_level is None:
                stderr_handler.setLevel("WARNING")
            else:
                stderr_handler.setLevel(args.log_level)
        else:
            log_file_handler = logging.FileHandler(filename=args.log, mode=mode)
            if args.log_level is None:
                log_file_handler.setLevel("INFO")
            else:
                log_file_handler.setLevel(args.log_level)
            log_file_handler.setFormatter(log_formatter)
            logger.addHandler(log_file_handler)
    return logger


def configure_subprocess_logger(queue):
    global logger
    logger = logging.getLogger(__name__)
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    logger.addHandler(QueueHandler(queue))
    logger.setLevel(logging.DEBUG)
    return logger


#class ForkedPdb(pdb.Pdb):
#    """A Pdb subclass that may be used
#    from a forked multiprocessing child
#
#    https://stackoverflow.com/questions/4716533/how-to-attach-debugger-to-a-python-subproccess/23654936#23654936
#    """
#    def interaction(self, *args, **kwargs):
#        _stdin = sys.stdin
#        try:
#            sys.stdin = open('/dev/stdin')
#            pdb.Pdb.interaction(self, *args, **kwargs)
#        finally:
#            sys.stdin = _stdin


def collapse_to_iupac(seqs):
    """Combine sequences into a consensus using IUPAC ambiguity codes
    
    Parameters:
    -----------
    seqs : list of str
        The sequences to combine.
    
    Returns
    -------
    str
        The combined sequence
    """
    seq_lens = [len(x) for x in seqs]
    max_len = max(seq_lens)
    if len(set(seq_lens)) != 1:  # TODO: replace with alignment
        return '-' * max_len
    output = []
    for i in range(max_len):
        column = {s[i] for s in seqs}
        if "*" in column or "N" in column or UNKNOWN_CHAR in column:
            output.append('N')
        else:
            output.append(iupac_key[tuple(sorted(column))])
    return "".join(output)


class GroupedRegion:

    def __init__(self, variants, group, reference, upstream=None, downstream=None):
        """
        Parameters
        ----------
        variants : iterator returning GroupedVariant
            Consecutive variants representing a region.
        group : str
            The group ID for the group to use.
        reference : str
            The reference sequence use to infer the variants
        upstream : deque of GroupedVariant
            The variants upstream (larger reference position) of the region of the interest
        downstream : deque of GroupedVariant
            The variants downstream (smaller reference position) of the region of the interest
       """
        if downstream is None:
            downstream = deque()
        if upstream is None:
            upstream = deque()
        self.variants = deque(variants)
        self.group = group
        self.reference = reference
        self.upstream = upstream
        self.downstream = downstream
        self.type = 'Undetermined'

    @classmethod
    def sliding_window(cls, variants, groups, reference, span, flank=1000):
        """Generate GroupedRegion objects as a sliding window along variants

        Parameters
        ----------
        variants : iterator returning GroupedVariant
            The variants to generate GroupedRegion objects from.
        groups : list of str
            The groups to return sliding windows for
        reference : str
            The file path to the reference sequence use to infer the variants
        span : int
            The number of nucleotides that the variants should span within a group
        flank : int
            The number of variants to keep in the upstream and downstream queues
        """

        def increment(region, flank):
            # Move one variant to spacer and add next to upstream queue
            region.variants.append(region.upstream.popleft())
            # Move from spacer to downstream queue once spacer is too large
            while region.region_length() > span:
                region.downstream.appendleft(region.variants.popleft())
            # Clear downstream queue once it is too large
            while len(region.downstream) > flank:
                region.downstream.pop()

        # Initialize sliding windows
        windows = {}
        for group in groups:
            windows[group] = cls(variants=[], group=group, reference=reference)
        # Read through variants and put into windows
        for index, variant in enumerate(variants):
            for group in groups:
                windows[group].upstream.append(variant)
                if index + 1 >= flank:  # Once the upstream is full
                    increment(windows[group], flank=flank)
                    if len(windows[group].variants) > 0:  # Can be 0 when a single var is larger than the span
                        yield cls(variants=windows[group].variants, group=group, reference=reference,
                                  upstream=windows[group].upstream, downstream=windows[group].downstream)
        # Clear the upstream window
        for index in range(len(next(iter(windows.values())).upstream)):
            for group in groups:
                increment(windows[group], flank=flank)
                if len(windows[group].variants) > 0:  # Can be 0 when a single var is larger than the span
                    yield cls(variants=windows[group].variants, group=group, reference=reference,
                              upstream=windows[group].upstream, downstream=windows[group].downstream)


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

    def region_length(self):
        """Get the maximum sequence length spanned by variants in a subset
        
        Return
        ------
        int/None
            The number of nucleotides spanned by the variants. If any
            of the samples have variable length alleles, then a single
            length cannot be determined and None is returned.
        """
        # Return 0 if there are no variants
        if len(self.variants) == 0:
            return 0

        # Get length of the reference spanned by these variants
        starts = [x.variant.pos - 1 for x in self.variants]
        ends = [s + x.variant.rlen - 1 for s, x in zip(starts, self.variants)]
        out = max(ends) - min(starts) + 1

        # Apply length changes for each variant to adjust length
        for var in self.variants:
            allele_lens = var.allele_lens(self.group)
            if len(allele_lens) == 0:
                continue
            max_allele_len = max(allele_lens.values())
            ref_len = len(var.variant.ref)
            out += max_allele_len - ref_len

        return out

    @staticmethod
    def _subset_sequence(variants, subset,
                         reference=None,
                         counts=True,
                         hetero=True,
                         unknown=True,
                         min_reads=0,
                         min_geno_qual=0,
                         min_samples=0,
                         consensus=False,
                         conserved=False,
                         wrap='<>',
                         sep=','):
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
        sep The character used to separate multiple alleles.
        
        Returns
        -------
        list of str
            A representation of a sequence, with one element per
            variant or single base from the reference sequence.
        """
        # Create rendered text for variant information
        rendered = []
        for variant in variants:
            allele_counts = GroupedVariant._count_genotypes(variant,
                                                            subset=subset,
                                                            hetero=hetero,
                                                            unknown=unknown,
                                                            min_reads=min_reads,
                                                            min_geno_qual=
                                                            min_geno_qual)
            sample_count = GroupedVariant._subset_sample_counts(variant, subset,
                                                                min_reads=min_reads,
                                                                min_geno_qual=
                                                                min_geno_qual)
            if len(allele_counts) == 0:
                allele_counts = {'?': 0}
            if conserved and (len(allele_counts) != 1 or sample_count <
                              min_samples):
                text = None
            else:
                if consensus:
                    text = collapse_to_iupac(allele_counts.keys())
                else:
                    if counts:
                        text_parts = [a + str(c) for a, c in
                                      allele_counts.items()]
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
            var_pos = [var.pos - 1 for var in variants]
            var_chrom = {var.chrom for var in variants}
            if len(var_chrom) > 1:
                raise ValueError('Variants cannot span multiple chromosomes')
            chrom = list(var_chrom)[0]
            ref_start = min(var_pos)
            ref_end = max(var_pos)
            ref_seq = list(reference[chrom].seq[ref_start:ref_end])
            # Replace reference with each variant
            for var, rend in zip(variants, rendered):
                replace_start = var.pos - 1 - ref_start
                replace_end = replace_start + len(var.ref)
                ref_seq = ref_seq[:replace_start] + [rend] + \
                          ref_seq[replace_end:]
            return ref_seq

    def sequence(self, reference, start, end, group=None, annotate=False):
        """Infer the sequence for the current group based on the reference sequence and variants
        """

        # Find variants within target region
        all_vars = self.downstream + self.variants + self.upstream
        var_starts = [x.variant.pos - 1 for x in all_vars]
        var_ends = [x.variant.pos + x.variant.rlen - 2 for x in all_vars]
        vars_in_range = [v for v, vs, ve in zip(all_vars, var_starts, var_ends) if start <= ve <= end or start <= vs <= end]

        # Check that variants are all on the same chormosome
        if len({x.variant.chrom for x in vars_in_range}) > 1:
            raise ValueError('Variants cannot span multiple chromosomes')
        chrom = self.variants[-1].variant.chrom

        # If there are no variants in the region, return the reference
        if len(vars_in_range) == 0:
            return list(reference[chrom].seq[start:end+1].lower())

        # Insert consensus for each variant in reference sequence
        # for var in vars_in_range:
        #     replace_start = var.variant.pos - 1 - start
        #     replace_end = replace_start + len(var.variant.ref)
        #     if group is None:
        #         consensus = var.variant.ref
        #     else:
        #         alleles = var.allele_counts[group].keys()
        #         if len(alleles) == 0:  # If no data for this position, use N of equal length to reference
        #             consensus = "N" * var.variant.rlen
        #         else:
        #             consensus = collapse_to_iupac(alleles)
        #     if var_upper:
        #         consensus = consensus.upper()
        #     ref_seq = ref_seq[:replace_start] + [consensus] + ref_seq[replace_end:]

        # Sort variants
        var_starts = [v.variant.pos - 1 for v in vars_in_range]
        var_ends = [s + v.variant.rlen - 1 for v, s in zip(vars_in_range, var_starts)]
        var_lens = [e - s + 1 for s, e in zip(var_starts, var_ends)]
        # var_diffs = [l - v.max_allele_len(group) for v, l in zip(vars_in_range, var_lens)]
        vars_in_range = [x for _, x in sorted(zip(var_ends, vars_in_range), key=lambda pair: pair[0])]

        # Adjust start and end positions of sequence to span all variants
        seq_ref_start = min(var_starts + [start])
        seq_ref_end = max(var_ends + [end])

        # Apply variants to reference sequence in reverse order
        out_seq = list(reference[chrom].seq[seq_ref_start:seq_ref_end+1].lower())
        for var in reversed(vars_in_range):
            replace_start = var.variant.pos - 1 - seq_ref_start
            replace_end = replace_start + len(var.variant.ref)
            is_diag_site = any([x is not None for x in var.diagnostic.values()])
            if group is None:
                replacement = var.variant.ref
            else:
                is_diag_for_group = var.diagnostic[group] is not None
                alleles = var.allele_counts[group]
                if len(alleles) == 0:  # If no data for this position, use N of equal length to reference
                    consensus = "N" * var.variant.rlen
                else:
                    consensus = collapse_to_iupac(alleles.keys())
                if annotate:
                    if is_diag_site:
                        replacement = ";".join([k + str(v) for k, v in alleles.items()])
                        if is_diag_for_group:
                            replacement = '<' + replacement + '>'
                    else:
                        replacement = consensus.upper()
                else:
                    if is_diag_for_group:
                        replacement = consensus.upper()
                    else:
                        replacement = consensus.lower()
            if annotate:
                out_seq = out_seq[:replace_start] + [replacement] + out_seq[replace_end:]
            else:
                out_seq = out_seq[:replace_start] + list(replacement) + out_seq[replace_end:]

        # Trim sequence to original start and end points
        if seq_ref_end > end:
            out_seq = out_seq[:len(out_seq) - (seq_ref_end - end)]
        if seq_ref_start < start:
            out_seq = out_seq[start - seq_ref_start:]

        return out_seq

    def conserved(self):
        """For each variant return alleles conserved in a given group

        Returns
        -------
        list of str/None
            Alleles for each variant when there is a conserved allele,
            otherwise None.
        """
        return [x.conserved[self.group] for x in self.variants]

    def diagnostic(self):
        """For each variant return alleles diagnostic for the group

        Returns
        -------
        list of str/None
            Alleles for each variant when there is a conserved allele,
            otherwise None.
        """
        return [x.diagnostic[self.group] for x in self.variants]

    def ref_pos_from_group_offset(self, ref_pos, offset):
        """Get ref index offset by an amount of group-specific sequence"""
        ref_diff_offset = 0
        for v in itertools.chain(reversed(self.downstream), self.variants, self.upstream):
            var_pos_diff = v.variant.pos - 1 - ref_pos
            var_group_offset = var_pos_diff + ref_diff_offset
            if var_group_offset >= offset:
                break
            if var_pos_diff >= 0:
                group_allele_len = v.max_allele_len(self.group)
                ref_allele_len = len(v.variant.ref)
                ref_diff_offset += group_allele_len - ref_allele_len
        return ref_pos + offset - ref_diff_offset


def _parse_reference(ref_path):
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


def parse_primer3_settings(file_path):
    """
    Reads primer3 BoulderIO format, assuming only global settings (starts with PRIMER_),
    and returns a dict in the format used by primer3-py.
    """
    def to_number_if_can(x):
        try:
            if int(float(x)) == float(x) and '.' not in x:
                return int(x)
            else:
                return float(x)
        except ValueError:
            return x

    with open(file_path) as handle:
        options = dict([tuple(l.strip().split('=')) for l in handle.readlines()])
        for opt, val in options.items():
            if ' ' in val or ';' in val:
                val = re.split('[ ;]+', val)
                val = [to_number_if_can(v) for v in val]
                if ',' in val or '-' in val[0]:
                    val = [[to_number_if_can(x) for x in re.split('[,-]+', v)] for v in val]
            elif ',' in val or '-' in val:
                val = re.split('[,\-]+', val)
                val = [to_number_if_can(v) for v in val]
            else:
                val = to_number_if_can(val)
            options[opt] = val
    return options


def run_primer3(template, target_start, target_len,
                options=None,
                tm=(53, 68),
                gc=(40, 70),
                amp_size=(80, 300),
                primer_size=(25, 35),
                max_sec_tm=40,
                gc_clamp=1,
                max_end_gc=4):
    if options is None:
        global_options = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,  # 1 == True
            # 'PRIMER_PICK_INTERNAL_OLIGO': 1,  # 1 == True
            'PRIMER_PICK_RIGHT_PRIMER': 1,  # 1 == True
            'PRIMER_LIBERAL_BASE': 1,  # 1 == True
            'PRIMER_OPT_SIZE': mean(primer_size),
            'PRIMER_MIN_SIZE': primer_size[0],
            'PRIMER_MAX_SIZE': primer_size[1],
            # 'PRIMER_INTERNAL_MAX_SIZE': len(crrna_seq),
            'PRIMER_OPT_TM': mean(tm),
            'PRIMER_MIN_TM': tm[0],
            'PRIMER_MAX_TM': tm[1],
            'PRIMER_MIN_GC': gc[0],
            'PRIMER_MAX_GC': gc[1],
            'PRIMER_MAX_POLY_X': 4,  # The maximum allowable length of a mononucleotide repeat
            'PRIMER_MAX_NS_ACCEPTED': 0,  # The maximum number of Ns
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,  # 1 == True
            'PRIMER_MAX_SELF_ANY_TH': max_sec_tm,  # describes the tendency of a primer to bind to itself
            'PRIMER_MAX_SELF_END_TH': max_sec_tm,  # 3' binding to itself (primer dimers)
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': max_sec_tm,  # the tendency of the left primer to bind to the right primer
            'PRIMER_PAIR_MAX_COMPL_END_TH': max_sec_tm,  # primer heterodimers
            'PRIMER_MAX_HAIRPIN_TH': max_sec_tm,  # hairpins
            'PRIMER_PRODUCT_SIZE_RANGE': [amp_size],
            'PRIMER_GC_CLAMP': gc_clamp,
            'PRIMER_MAX_END_GC': max_end_gc,
        }
    else:
        global_options = parse_primer3_settings(options)

    p3_output = primer3.bindings.design_primers(
        {
            'SEQUENCE_TEMPLATE': "".join(template),
            # 'SEQUENCE_INTERNAL_OLIGO': "".join(crrna_seq),
            'SEQUENCE_TARGET': [target_start, target_len]
        },
        global_options
    )
    return p3_output


def consv_border_n(group, border_var, nearby_vars, max_offset):
    """Find the maximum size of adjacent conserved reference sequence in terms of reference and group sequence"""
    # Initialize cumulative length differences between reference and variant alleles
    if len(nearby_vars) > 0 and border_var.variant.pos < nearby_vars[0].variant.pos:
        ref_diff_offset = border_var.max_allele_len(group) - len(border_var.variant.ref)
    else:
        ref_diff_offset = 0
    # Initialize distance between reference border var and the start/end of a variant
    ref_diff = 0
    # Search though variants until a non-conserved variant is found
    for nearby_var in nearby_vars:
        # Get info about reference and group
        group_len = nearby_var.max_allele_len(group)
        ref_len = len(nearby_var.variant.ref)
        ref_start = nearby_var.variant.pos
        ref_end = ref_start + ref_len - 1
        # Get distance between reference border var and the start/end of a variant
        if border_var.variant.pos <= ref_start:
            ref_diff = ref_start - border_var.variant.pos
        else:
            ref_diff = border_var.variant.pos - ref_end
        # Check if variant is further away than the maximum offset
        if ref_diff + ref_diff_offset >= max_offset:
            return {"ref": max_offset - ref_diff_offset, "group": max_offset}
        # If variant is not conserved, return pos right before variant
        if nearby_var.conserved[group] is None:
            return {"ref": ref_diff - 1,
                    "group": ref_diff + ref_diff_offset - 1}
        # Keep track of differences in allele length between the group and the reference
        ref_diff_offset += group_len - ref_len
    return {"ref": ref_diff - ref_diff_offset, "group": ref_diff}


def is_nearby_conserved(group, border_var, nearby_vars, max_offset):
    """Check if some number of bases near a variant are conserved """
    ref_diff_offset = 0  # cumulative differences between reference and variant alleles
    for nearby_var in nearby_vars:
        max_allele_len = nearby_var.max_allele_len(group)
        ref_diff_offset += max_allele_len - len(nearby_var.variant.ref)
        var_pos_diff = abs(border_var.variant.pos - nearby_var.variant.pos)
        if var_pos_diff + ref_diff_offset > max_offset:
            return True
        if nearby_var.conserved[group] is None:
            return False


class DiagosticRegion(GroupedRegion):
    """The class for reporting diagnostic regions."""

    def __init__(self,  variants, group, reference, upstream, downstream, p3, crrna_seq, downstream_seq, upstream_seq, temp_range, crrna_range, min_bases):
        super().__init__(variants, group, reference, upstream, downstream)
        self.p3 = p3
        self.downstream_seq = downstream_seq
        self.crrna_seq = crrna_seq
        self.upstream_seq = upstream_seq
        self.temp_range = temp_range
        self.crrna_range = crrna_range
        self.type = "Diagnostic"
        self.min_bases = min_bases

    @staticmethod
    def from_grouped_region(region, p3, crrna_seq, downstream_seq, upstream_seq, temp_range, crrna_range, min_bases):
        return DiagosticRegion(variants=region.variants,
                               group=region.group,
                               reference=region.reference,
                               upstream=region.upstream,
                               downstream=region.downstream,
                               p3=p3,
                               crrna_seq=crrna_seq,
                               downstream_seq=downstream_seq,
                               upstream_seq=upstream_seq,
                               temp_range=temp_range,
                               crrna_range=crrna_range,
                               min_bases=min_bases)

    def left_range(self):
        """Start/stop reference sequence indexes of the left primer."""
        start = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=self.p3['PRIMER_LEFT_0'][0])
        end = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=sum(self.p3['PRIMER_LEFT_0']) - 1)
        return [start, end]

    def right_range(self):
        """Start/stop reference sequence indexes of the right primer."""
        start = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0],
                                               offset=self.p3['PRIMER_RIGHT_0'][0] - self.p3['PRIMER_RIGHT_0'][1] + 1)
        end = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0],
                                             offset=self.p3['PRIMER_RIGHT_0'][0])
        return [start, end]

    def missing_samples(self):
        diag_vars = [var for var in self.variants if var.diagnostic[self.group] is not None]
        return {id for var in diag_vars for ids in var.missing_samp_ids.values() for id in ids}

    # def crrna_range(self):
    #     """Start/stop reference sequence indexes of the crRNA."""
    #     # start = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=self.p3['PRIMER_INTERNAL_0'][0])
    #     # end = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=sum(self.p3['PRIMER_INTERNAL_0']))
    #     # return [start, end]
    #     return self.crrna_range


def find_diag_region(variants,
                     groups,
                     reference=None,
                     nontarget=None,
                     primer3=False,
                     min_vars=1,
                     min_bases=1,
                     min_samp_prop=0.9,
                     min_samples=5,
                     min_reads=5,
                     min_geno_qual=30,
                     min_map_qual=40,
                     min_var_qual=10,
                     min_freq=0.1,
                     crrna_len=28,
                     tm = (53, 68),
                     gc = (40, 70),
                     amp_size=(80, 300),
                     primer_size=(25,35),
                     max_sec_tm = 40,
                     gc_clamp=1,
                     max_end_gc=4,
                     var_location=(4,16),
                     force=False):
    """Find regions with diagnostic variants
    
    Return information on regions that contain variants diagnostic for a subset
    of samples, optionally surrounded by conserved regions where primers can be
    designed.
    
    Parameters
    ----------
    variants : an iterable returning variants or str
        A series of consecutive variants in which to identify diagnostic
        clusters of variants.
    groups : dict of list/tuple of str
        The sample IDs for each group to find diagnostic clusters for,
        named by group.
    nontarget : list of str, optional
        The sample IDs that should be distinct from those in `groups`, but
        are otherwise not of interest. `min_samples` does not apply to these.
    reference :
        The path to a FASTA file. Required if Primer3 is to be used.
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
    min_samp_prop : float, optional
        The minimum proportion of samples that must contain a diagnositc variant in each group
    min_samples : int, optional
        The minimum number of samples that must represent each group in
        `groups` for a given variants. Samples must pass the `min_reads` 
        and `min_geno_qual` filters to count towards this minimum.
    min_reads : int, optional
        The minimum number of reads a sample have at the location of a
        given variant.
    min_geno_qual : int, optional
        The minimum genotype quality score (phred scale). This corresponds
        to the per-sample output in the VCF encoded as "GQ".
    min_map_qual : int, optional
        The minimum root square mean mapping quality of reads supporting a
        variant (phred scale). This corresponds to the per-variant
        output in the VCF encoded as "MQ".
    min_var_qual : int, optional
        Phred-scaled quality score for the assertion that a variant exists.
        This corresponds to the per-variant output in the VCF encoded as the "QUAL" column.
    min_freq : float, optional
        The minimum proportion of reads an allele must have to be considered real.
        This is meant to counter sequencing errors.
        If not set, the genotypes called in the VCF are used without considering read depth.
    crrna_len: int, optional
        The length of the spacer sequence that the crRNA will bind to and
        the primers will amplify.
    offset_left : int, optional
        The minimum number of bases on the 3' end of the spacer before
        the first diagnostic SNP.
    offset_right : int, optional
        The number of bases on the 5' end after the last diagnostic SNP.
        
    Yields
    ------
    object?
        information on regions that contain variants diagnostic for a
        subset of samples. TBD.
    """
    offset_left = var_location[0] - 1
    offset_right = crrna_len - var_location[1]
    window_width = crrna_len - offset_right - offset_left

    vcf_reader = GroupedVariant.from_vcf(variants, groups,
                                         min_samp_prop=min_samp_prop,
                                         min_samples=min_samples,
                                         min_reads=min_reads,
                                         min_geno_qual=min_geno_qual,
                                         min_freq=min_freq,
                                         min_map_qual=min_map_qual,
                                         min_var_qual=min_var_qual,
                                         force=force)
    windower = GroupedRegion.sliding_window(vcf_reader,
                                            groups=groups.keys(),
                                            reference=reference,
                                            span=window_width,
                                            flank=amp_size[1])
    for region in windower:

        # Are there enough diagnostic variants?
        is_diag = [x is not None for x in region.diagnostic()]
        n_diag_var = sum(is_diag)
        if n_diag_var < min_vars:
            region.type = 'Undiagnostic'
            yield region
            continue

        # Are all the variants in the crRNA conserved?
        if any([x is None for x in region.conserved()]):
            region.type = 'Unconserved'
            yield region
            continue

        # If there is only one diagnostic variants, is it on the right?
        # NOTE: this might miss some valid regions when a conservered variant is a few bp more to the right
        if n_diag_var == 1 and is_diag[-1] is False:
            region.type = 'Misplaced'
            yield region
            continue

        # # Is there conserved sequence adjacent to the variants to fit the crRNA?
        # region_len = region.region_length()
        # overhang_left = spacer_len - offset_right - region_len
        # if not is_nearby_conserved(region.group, region.variants[-1], region.upstream, offset_right + 1):
        #      region.type = 'Unconserved'
        #      yield region
        # if not is_nearby_conserved(region.group, region.variants[0], region.downstream, overhang_left + 1):
        #      region.type = 'Unconserved'
        #      yield region

        # Is there conserved adjacent sequence for the crRNA?
        overhang_left = crrna_len - region.region_length() - offset_right
        overhang_right = offset_right
        overhang_len_up = consv_border_n(group=region.group,
                                         border_var=region.variants[-1],
                                         nearby_vars=region.upstream,
                                         max_offset=overhang_right)
        overhang_len_dn = consv_border_n(group=region.group,
                                         border_var=region.variants[0],
                                         nearby_vars=region.downstream,
                                         max_offset=overhang_left)
        if overhang_len_up['group'] < offset_right or overhang_len_dn['group'] < overhang_left:
            region.type = 'Unconserved'
            yield region
            continue

        # Is there enough adjacent conserved regions to design primers?
        consv_len_up = consv_border_n(group=region.group,
                                      border_var=region.variants[-1],
                                      nearby_vars=region.upstream,
                                      max_offset=amp_size[1])
        consv_len_dn = consv_border_n(group=region.group,
                                      border_var=region.variants[0],
                                      nearby_vars=region.downstream,
                                      max_offset=amp_size[1])
        if consv_len_up["group"] - overhang_len_up['group'] < 30:  #TODO base on primer3 primer size
            region.type = 'Unconserved'
            yield region
            continue
        if consv_len_dn["group"] - overhang_len_dn['group'] < 30:  #TODO base on primer3 primer size parameters
            region.type = 'Unconserved'
            yield region
            continue


        # Are there enough diagnostic bases?
        start_crrna_ref = region.variants[0].variant.pos - 1 - overhang_len_dn['ref']
        end_crrna_ref = region.variants[-1].variant.pos - 1 + overhang_len_up['ref']
        crrna_seq = region.sequence(reference=reference, start=start_crrna_ref, end=end_crrna_ref, group=region.group)
        crrna_nontaget_seqs = [region.sequence(reference=reference, start=start_crrna_ref, end=end_crrna_ref, group=g) for g in groups if g is not region.group]
        edit_dists = [distance.edit_distance("".join(crrna_seq), "".join(seq)) for seq in crrna_nontaget_seqs]
        min_dist = min(edit_dists)
        if min_dist < min_bases:
            region.type = 'Undiagnostic'
            yield region
            continue

        start_tmp_ref = region.variants[0].variant.pos - 1 - consv_len_dn['ref']
        end_tmp_ref = region.variants[-1].variant.pos - 1 + consv_len_up['ref']

        # Are there no overlapping variants TODO: fix so this is not needed
        #all_vars = region.downstream + region.variants + region.upstream
        #var_starts = [x.variant.pos - 1 for x in all_vars]
        #var_ends = [x.variant.pos + x.variant.rlen - 2 for x in all_vars]
        #vars_in_range = [v for v, vs, ve in zip(all_vars, var_starts, var_ends) if start_tmp_ref <= ve <= end_tmp_ref or start_tmp_ref <= vs <= end_tmp_ref]
        #overlapping = False
        #for v1, v2 in zip(list(vars_in_range)[:-1], list(vars_in_range)[1:]):
        #    if v1.variant.stop > v2.variant.start:
        #        region.type = 'Overlapping'
        #        overlapping = True
        #        break
        #if overlapping:
        #    yield region
        #    continue

        # Run primer3 on group-specific template
        downstream_seq = region.sequence(reference=reference, start=start_tmp_ref, end=start_crrna_ref-1, group=region.group)
        upstream_seq = region.sequence(reference=reference, start=end_crrna_ref+1, end=end_tmp_ref, group=region.group)

        template_seq = downstream_seq + crrna_seq + upstream_seq
        start_crrna_tmp = len(downstream_seq)
        end_crrna_tmp = start_crrna_tmp + len(crrna_seq) - 1

        p3_out = run_primer3(template_seq, target_start=start_crrna_tmp, target_len=len(crrna_seq),
                             tm=tm,
                             gc=gc,
                             amp_size=amp_size,
                             primer_size=primer_size,
                             max_sec_tm=max_sec_tm,
                             gc_clamp=gc_clamp,
                             max_end_gc=max_end_gc)

        if p3_out['PRIMER_PAIR_NUM_RETURNED'] == 0:
            region.type = 'No primers'
            yield region
            continue

        # Return data for the diagnostic region
        region.type = 'Diagnostic'
        output = DiagosticRegion.from_grouped_region(region, p3=p3_out, crrna_seq=crrna_seq, downstream_seq=downstream_seq, upstream_seq=upstream_seq,
                                                     temp_range=[start_tmp_ref, end_tmp_ref], crrna_range=[start_crrna_ref, end_crrna_ref], min_bases=min_dist)

        yield output


def parse_command_line_args():
    # Parse command line arguments
    logger.debug('Parsing command line arguments')
    parser = argparse.ArgumentParser(
        description='Find regions where there are conserved variants for each group that are not found in other groups.')
    parser.add_argument('metadata', type=str, metavar='METADATA',
                        help='A CSV file containing data with one row per sample. Two columns are required: one which contains the same sample IDs used in the VCF file, and one which identifies which group each sample belongs to. See --sample_col and --group_col for default values and how to specify custom column names.')
    parser.add_argument('reference', type=str, metavar='REFERENCE',
                        help='The reference file used to make the VCF input. If supplied, the sequence of the region containing each conserved variant is returned with the output.')
    parser.add_argument('--vcf', type=str, default="-", metavar='PATH',
                        help='A VCF file containing variant data for the samples grouped by the `metadata` file. If not supplied, VCF data will be read from stanard input (stdin) using a single core. (default: read from stdin)')
    parser.add_argument('--sample_col', type=str, default="sample_id", metavar='TEXT',
                        help='The names of column in the metadata containing sample IDs. (default: %(default)s)')
    parser.add_argument('--group_col', type=str, default="group", metavar='TEXT',
                        help='The names of column in the metadata containing group IDs. Samples with the same groups ID will represent that group. (default: %(default)s)')
    parser.add_argument('--index', type=str, metavar='PATH',
                        help='The path to an tabix index file for the VCF file. If not supplied, a file with the same name as the VCF file with .tbi/.csi appended will be searched for. If that is not found, an index file will be created in the same directory as the VCF file.')
    parser.add_argument('--groups', type=str, nargs="+", metavar='TEXT',
                        help='One or more groups that are to be distinguished by variants. These should match the values of the column specified by --group_col in the metadata file. (default: use all groups)')
    parser.add_argument('--out_csv', type=str, metavar='PATH',
                        help='The output file to create. If not supplied, results will be printed to the screen (standard out). (default: print to stdout)')
    parser.add_argument('--out_align', type=str, metavar='PATH',
                        help='A file path to print human-readable alignments of diagnostic regions. (default: do not output)')
    parser.add_argument('--chroms', type=str, nargs="+", metavar='TEXT',
                        help='One or more chromosomes (reference fasta headers) to restrict the search to. (default: use all chromosomes)')
    parser.add_argument('--pos', type=int, nargs=2, metavar='INT', default=None,
                        help='The range of indexes (1-based) of positions in the reference sequences to search. (default: search whole chromosome)')
    parser.add_argument('--min_samples', type=int, default=3, metavar='INT',
                        help='The number of samples with data passing quality filters each group must have for a given variant. (default: %(default)s)')
    parser.add_argument('--min_samp_prop', type=float, default=0.9, metavar='PROP',
                        help='The minimum proportion of samples that must contain a diagnositc variant in each group. (default: %(default)s)')
    parser.add_argument('--min_reads', type=int, default=10, metavar='INT',
                        help='The number of reads a variant must be represented by in a given sample for the data of that sample to be considered. This corresponds to the per-sample GATK output in the VCF encoded as "DP". (default: %(default)s)')
    parser.add_argument('--min_geno_qual', type=int, default=40, metavar='INT',
                        help='The minimum genotype quality score (phred scale). This corresponds to the per-sample GATK output in the VCF encoded as "GQ". (default: %(default)s)')
    parser.add_argument('--min_var_qual', type=int, default=10, metavar='INT',
                        help='Phred-scaled quality score for the assertion that a variant exists. This corresponds to the per-variant output in the VCF encoded as the "QUAL" column. (default: %(default)s)')
    parser.add_argument('--min_freq', type=float, default=0.1, metavar='PROP',
                        help='The minimum proportion of reads an allele must have to be considered real. This is meant to counter sequencing errors. (default: %(default)s)')
    parser.add_argument('--min_map_qual', type=int, default=40, metavar='INT',
                        help='The minimum root square mean mapping quality of reads supporting a variant (phred scale). This corresponds to the per-variant output in the VCF encoded as "MQ" in the "INFO" column. (default: %(default)s)')
    parser.add_argument('--min_bases', type=int, default=1, metavar='INT',
                        help='The minimum number of bases distinguishing the target group from other groups. (default: %(default)s)')
    parser.add_argument('--cores', type=int, default=1, metavar='INT',
                        help='The number of processors to use for parallel processing. (default: %(default)s)')
    parser.add_argument('--log', type=str, metavar='PATH',
                        help='The location to save a log file containing information, warnings, and errors. (default: print to screen via stderr)')
    parser.add_argument('--log_level', type=str, choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='The minimum importance level of messages to print to the log. (default: INFO if the log is being saved to a file with --log, otherwise WARNING)')
    parser.add_argument('--var_location', type=int, nargs=2, metavar='INT', default=[6, 14],
                        help="The range of acceptable locations for diagnostic variants in the crRNA, measured from the 5' end. (default: %(default)s)")
    parser.add_argument('--crrna_len', type=int, default=28, metavar='INT',
                        help='Then length of the crRNA. (default: %(default)s)')
    parser.add_argument('--tm', type=int, nargs=2, metavar='INT', default=[53, 68],
                        help='The minimum and maximum melting temperature when searching for primers. (default: %(default)s)')
    parser.add_argument('--gc', type=int, nargs=2, metavar='INT', default=[40, 70],
                        help='The minimum and maximum GC percentage when searching for primers. (default: %(default)s)')
    parser.add_argument('--amp_size', type=int, nargs=2, metavar='INT', default=[70, 150],
                        help='The minimum and maximum size of the amplicon when searching for primers. (default: %(default)s)')
    parser.add_argument('--primer_size', type=int, nargs=2, metavar='INT', default=[25, 35],
                        help='The minimum and maximum size of the primers. (default: %(default)s)')
    parser.add_argument('--max_sec_tm', type=int, default=40, metavar='INT',
                        help='The maximum melting temperature of any secondary structures when searching for primers, including hetero/homo dimers and hairpins. (default: %(default)s)')
    parser.add_argument('--gc_clamp', type=int, default=1, metavar='INT',
                        help="Require the specified number of consecutive Gs and Cs at the 3' end of both the left and right primer. (default: %(default)s)")
    parser.add_argument('--max_end_gc', type=int, default=4, metavar='INT',
                        help="The maximum number of Gs or Cs allowed in the last five 3' bases of a left or right primer. (default: %(default)s)")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Try to continue executing when the program would typically quit.")
    args = parser.parse_args()

    return args


def read_vcf_contigs(path, reference, index=None, chunk_size=100000, flank_size=1000, contig_subset=None, pos_subset=None):
    """Supply iterators for contigs/scaffolds/chormosomes in a VCF file.

    Parameters
    ----------
    path : str
        The path to the VCF file
    index : str, optional
        The path to the index file for the VCF
    """
    # If reading from stdin, return None
    if path == "-":
        return [None]
    # Get path to index file or create it
    if index is None:
        tbi_path = path + '.tbi'
        csi_path = path + '.csi'
        if os.path.isfile(tbi_path):
            index = tbi_path
        elif os.path.isfile(csi_path):
            index = csi_path
        else:
            logger.info(f'Creating index file  "{tbi_path}"')
            pysam.tabix_index(path, preset='vcf', keep_original=True, force=True)
            index = path + '.gz.tbi'
            path = path + '.gz'
    # Reduce chunk size if too big for custom position range
    if pos_subset is not None:
        pos_length = max(pos_subset) - min(pos_subset) + 1
        if pos_length < chunk_size:
            chunk_size = pos_length
    # Split contigs into chunks for processing
    index_handle = pysam.TabixFile(filename=path, index=index)
    output = []
    for contig in index_handle.contigs:
        if contig_subset is not None and contig not in contig_subset:
            continue
        if pos_subset is None:
            search_start = 0
            search_end = len(reference[contig])
        else:
            search_start = min(pos_subset) - 1
            search_end = max(pos_subset) - 1
        for start in range(search_start, search_end, chunk_size):
            end = start + chunk_size + flank_size
            if start > flank_size:
                start -= flank_size
            output.append({'contig': contig, 'start': start, 'end': end})

    return output


def _format_p3_output(p3_out):
    """Reformat data for best primer pair for CSV output"""
    return {primer3_col_key[n]: p3_out[n] for n in primer3_col_names}


def _render_header(stream, out_sep=','):
    print("group", "chrom", "n_diag", "fwd_start", "fwd_end", "rev_start",
          "rev_end", "crrna_start", "crrna_end", "seq_start", "seq_end",
          *[primer3_col_key[n] for n in primer3_col_names],
          sep=out_sep, file=stream)


@contextmanager
def stream_writer(file_path=None, default_stream=sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "w")
    yield file_handle
    if file_path is not None:
        file_handle.close()


def _format_for_csv(region, reference, groups):
    fwd_range = region.left_range()
    rev_range = region.right_range()
    crrna_range = region.crrna_range
    temp_range = region.temp_range

    group = region.group
    chrom = region.variants[0].variant.chrom
    # n_diag = sum([x is not None for x in region.diagnostic()])
    n_diag = region.min_bases

    fwd_start = fwd_range[0] + 1  # The +1 goes from 0-based in 1-based indexing
    fwd_end = fwd_range[1] + 1
    rev_start = rev_range[0] + 1
    rev_end = rev_range[1] + 1
    crrna_start = crrna_range[0] + 1
    crrna_end = crrna_range[1] + 1
    seq_start = region.temp_range[0] + 1
    seq_end = region.temp_range[1] + 1

    # Make unannotated sequences of target group
    def format_seq(start, end):
        out = region.sequence(start=start, end=end,
                              reference=reference, group=region.group, annotate=False)
        return "".join(out)

    seq_adj_left = format_seq(start=temp_range[0], end=fwd_range[0] - 1)
    seq_primer_left = format_seq(start=fwd_range[0], end=fwd_range[1])
    seq_amp_left = format_seq(start=fwd_range[1] + 1, end=crrna_range[0] - 1)
    seq_crrna = format_seq(start=crrna_range[0], end=crrna_range[1])
    seq_amp_right = format_seq(start=crrna_range[1] + 1, end=rev_range[0] - 1)
    seq_primer_right = format_seq(start=rev_range[0], end=rev_range[1])
    seq_adj_right = format_seq(start=rev_range[1] + 1, end=temp_range[1])

    missing = region.missing_samples()
    missing_samp_ids = ';'.join(missing)
    missing_count = len(missing)


    # Modify these values to change columns and their names:
    output = {
        "region_id": f'{chrom}:{fwd_start}-{rev_end}',
        "group": group,
        "chrom": chrom,
        "n_diag": n_diag,
        "n_missing": missing_count,
        "reg_from": seq_start,
        "reg_to": seq_end,
        "diag_from": crrna_start,
        "diag_to": crrna_end,
        "fwd_from": fwd_start,
        "fwd_to": fwd_end,
        "rev_from": rev_start,
        "rev_to": rev_end,
        "missing_samp_ids": missing_samp_ids,
        "seq_adj_left": seq_adj_left,
        "seq_primer_fwd": seq_primer_left,
        "seq_inter_left": seq_amp_left,
        "seq_diag": seq_crrna,
        "seq_inter_right": seq_amp_right,
        "seq_primer_rev": seq_primer_right,
        "seq_adj_right": seq_adj_right
    }
    output.update(_format_p3_output(region.p3))
    # output.update(group_seqs)

    return output


def _print_alignment(region, reference, groups):

    fwd_range = region.left_range()
    rev_range = region.right_range()
    crrna_range = region.crrna_range

    # Format sequences for each group
    def format_seq(group, start, end):
        out = region.sequence(start=start, end=end,
                              reference=reference, group=group, annotate=True)
        return out
    group_seqs = {g: format_seq(g, start=fwd_range[0], end=rev_range[1]) for g in groups}
    ref_seq = format_seq(None, start=fwd_range[0], end=rev_range[1])

    # Format primer and crrna sequences
    def format_oligo(start, end):
        out = region.sequence(start=start, end=end,
                              reference=reference, group=region.group, annotate=False)
        return "".join(out)
    seq_primer_left = format_oligo(start=fwd_range[0], end=fwd_range[1])
    seq_primer_right = format_oligo(start=rev_range[0], end=rev_range[1])
    seq_crrna = format_oligo(start=crrna_range[0], end=crrna_range[1])
    # oligos = {
    #     seq_primer_left: 0,
    #     seq_crrna: crrna_range[0] - fwd_range[0],
    #     seq_primer_right: rev_range[0] - fwd_range[0]
    # }
    oligos = [
        Annotation(name="Left primer", seq=seq_primer_left, start=0),
        Annotation(name="crRNA", seq=seq_crrna, start=crrna_range[0] - fwd_range[0]),
        Annotation(name="Right primer", seq=seq_primer_right, start=rev_range[0] - fwd_range[0])
    ]

    chrom = list(region.reference.keys())[0]
    start = fwd_range[0] + 1
    end = rev_range[1] + 1
    group = region.group
    output = [f"## {chrom}:{start}-{end} is diagnostic for {region.group}\n"]
    try:
        output += render_variant(seqs=group_seqs, ref=ref_seq, p3=region.p3, groups=groups, annots=oligos)
    except (IndexError, TypeError) as error:
        output += ["CANNOT PRINT ALIGNMENT WITH OVERLAPPING INDELS"]
        logger.info(f"Failed to print alignment of {chrom}:{start}-{end} diagnostic for {region.group} due to overlapping indels.")
    output += ['\n']

    return output


def report_diag_region(vcf_path, contig, groups, reference, args, **kwargs):
    if contig is None:  # If reading from stdin
        logger.info(f"Reading VCF data from standard input (stdin).")
        variants = pysam.VariantFile("-")
    else:
        logger.info(f"Starting scan of contig {contig['contig']}:{contig['start']}-{contig['end']}")
        vcf_handle = pysam.VariantFile(vcf_path)
        variants = vcf_handle.fetch(contig['contig'], start=contig['start'], end=contig['end'])
    stats = defaultdict(int)
    undiag_count = 0 # How many variants are not diagnositic and have therefore not been reported yet
    update_interval = 1000 # How often stats on non-diagnostic variants are reported
    for region in find_diag_region(variants, groups, reference, **kwargs):
        if failure_event is not None and failure_event.is_set():
            logger.critical("Error detected in other worker process. Ending this process too.")
            return None
        stats[region.type] += 1
        if region.type == "Diagnostic":
            output = _format_for_csv(region, reference, groups)
            if args.out_align is None:
                alignment = None
            else:
                alignment = _print_alignment(region, reference, groups)
            yield {'result': output, 'stats': stats, 'alignment': alignment}
            stats = defaultdict(int)
        else:
            undiag_count += 1
        if undiag_count >= update_interval:
            yield {'result': None, 'stats': stats, 'alignment': None}
            undiag_count = 0
            stats = defaultdict(int)
    return None




class ResultWriter:

    def __init__(self, output_stream, groups, align_path=None):
        self.result_header_printed = False
        self.stat_header_printed = False
        self.stats = defaultdict(int)
        self.output_stream = output_stream
        self.stat_names = ['Undiagnostic', 'Unconserved', 'No primers']
        self.variant_counts = {s: 0 for s in self.stat_names}
        self.groups = groups
        self.group_counts = {g: 0 for g in groups}
        self.align_path = align_path
        if align_path is not None:
            self.out_align = open(align_path, "w")

    def print_result(self, result):
        if not self.result_header_printed:
            print(*result.keys(), sep=',', file=self.output_stream, flush=True)
            self.result_header_printed = True
        print(*result.values(), sep=',', file=self.output_stream, flush=True)

    def print_stats_header(self):
        max_nchar = max([len(n) for n in self.stat_names + self.groups])
        header_parts = [n.ljust(max_nchar) for n in self.stat_names + self.groups]
        print('| '.join(header_parts), file=sys.stderr)

    def print_status(self, end_line=False):
        if not self.stat_header_printed:
            self.print_stats_header()
            self.stat_header_printed = True
        max_nchar = max([len(n) for n in self.stat_names + self.groups])
        var_info = [str(self.variant_counts[n]).ljust(max_nchar) for n in self.stat_names]
        group_info = [str(self.group_counts[n]).ljust(max_nchar) for n in self.groups]
        print('| '.join(var_info + group_info), file=sys.stderr, end='\n' if end_line else '\r')

    def update_stats(self, output):
        if output['result'] is not None:
            self.group_counts[output['result']['group']] += 1
        for stat, count in output['stats'].items():
            if stat in self.variant_counts:
                self.variant_counts[stat] += count

    def write_alignment(self, output):
        if self.align_path is not None:
            self.out_align.writelines([x + '\n' for x in output] + ['\n'])

    def write(self, output):
        if output['result'] is not None:
            self.print_result(output['result'])
            self.write_alignment(output['alignment'])
        self.update_stats(output)
        self.print_status()

    def finish(self):
        print('', file=sys.stderr) #move line down so stats are not covered


def mp_worker(queue, log_queue, args, contig, groups, reference, **kwargs):
    configure_subprocess_logger(log_queue)
    try:
        for result in report_diag_region(args.vcf, contig, groups, reference, args, **kwargs):
            queue.put(result)
    except BaseException as error:
        logger.exception(f"Error encountered while scanning contig {contig['contig']}:{contig['start']}-{contig['end']}:")
        failure_event.set()
        raise error
    return None


def mp_listener(result_queue, log_queue, args):
    '''listens for output on the queue and writes to a file. '''
    global logger
    logger = configure_global_logger(args, mode="a")
    with stream_writer(args.out_csv, sys.stdout) as output_stream:
        writer = ResultWriter(output_stream, args.groups, align_path=args.out_align)
        while True:
            # Write one result returned by workers
            try:
                output = result_queue.get(block=False, timeout=0.1)
                if output == "kill":
                    break
                writer.write(output)
            except queue.Empty:
                pass

            # Write all logs returned by workers
            while True:
                try:
                    logs = log_queue.get(block=False)
                    # logger = logging.getLogger(logs.name)
                    logger.handle(logs)
                except queue.Empty:
                    break
    writer.finish()
    total_vars = sum(writer.variant_counts.values()) + sum(writer.group_counts.values())
    logger.info("Total variants scanned: " + str(total_vars))


def mp_worker_init(event):
    global failure_event
    failure_event = event


def run_all():

    # Parse command line arguments
    args = parse_command_line_args()

    # Log record of parsed parameters
    global logger
    logger = configure_global_logger(args)
    lines = [f"    {k : <15}: {v}" for k, v in vars(args).items() if v is not None]
    logger.info("\n".join(["Parameters used:"] + lines))

    # Prepare input data
    reference = _parse_reference(args.reference)
    groups = _parse_group_data(args.metadata, groups=args.groups, sample_col=args.sample_col,
                               group_col=args.group_col, min_samples=args.min_samples)
    contigs = read_vcf_contigs(args.vcf, reference=reference, chunk_size=100000, flank_size=1000,
                               contig_subset=args.chroms, pos_subset=args.pos) #TODO base on amplicon size, need to add to args
    gz_path = args.vcf + '.gz'
    if os.path.isfile(gz_path):
        args.vcf = gz_path
    search_arg_names = ('min_samples', 'min_reads', 'min_geno_qual', 'min_map_qual', 'min_var_qual', 'min_freq',
                        'min_samp_prop', 'var_location', 'crrna_len', 'tm', 'gc', 'primer_size', 'amp_size',
                        'max_sec_tm', 'min_bases', 'gc_clamp', 'max_end_gc', 'force')
    search_args = {k: v for k, v in vars(args).items() if k in search_arg_names}

    if args.vcf != "-" and args.cores > 1:
        # Prepare for multiprocessing
        logger.debug('Preparing for multiprocessing')
        mp.set_start_method('spawn')
        manager = mp.Manager()
        failure_event = manager.Event()
        queue = manager.Queue()
        log_queue = manager.Queue()
        pool = mp.Pool(args.cores + 1, initializer=mp_worker_init, initargs=(failure_event,))
        watcher = pool.apply_async(mp_listener, args=(queue, log_queue, args))

        # Set main process logger to submit to the listener process instead
        configure_subprocess_logger(log_queue)

        # Spawn processes
        logger.debug('Spawning multiprocessing workers')
        jobs = []
        for contig in contigs:
            job = pool.apply_async(mp_worker,
                                   args=(queue, log_queue, args, contig, groups, reference),
                                   kwds=search_args)
            jobs.append(job)

        # End multiprocessing
        logger.debug('Waiting for multiprocessing workers to finish')
        for job in jobs:
            job.get()
        queue.put('kill')
        pool.close()
        pool.join()

        # Reset logger to use main process
        logger = configure_global_logger(args, mode='a')
    else:
        logger.debug('Running code on a single core')
        with stream_writer(args.out_csv, sys.stdout) as output_stream:
            writer = ResultWriter(output_stream, args.groups, align_path=args.out_align)
            for contig in contigs:
                for result in report_diag_region(args.vcf, contig, groups, reference, args,
                                                 **search_args):
                    writer.write(result)
            writer.finish()
        total_vars = sum(writer.variant_counts.values()) + sum(writer.group_counts.values())
        logger.info("Total variants scanned: " + str(total_vars))


def main():
    # import cProfile; cProfile.runctx('run_all()', globals=globals(), locals=locals(), sort='cumulative')
    run_all()


if __name__ == "__main__":
    # Test command:
    # python -m diagvar.find_diag_region test_data/test_metadata.csv test_data/unfilt_allscafs_n666.vcf.gz --groups NA1 NA2 EU1 EU2 --reference test_data/PR-102_v3.1_s0001.fasta --out test.csv
    main()

# if __name__ == "__main__":
#     vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
#     groups = _read_group_data('test_data/test_metadata.csv')
#     groups = {g: v for g, v in groups.items() if g in ["NA1", "NA2", "EU1", "EU2"]}
#     ref = 'test_data/PR-102_v3.1.fasta'
#     find_diag_region(vcf,
#                      groups=groups,
#                      reference=ref)

