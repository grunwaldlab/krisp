import sys
import argparse
import pandas
import gzip
import copy
import primer3
import pysam
import re
import itertools
from collections import deque, defaultdict
from Bio import SeqIO
from Bio import Seq
from Bio.Data import IUPACData
from contextlib import contextmanager
from collections import Counter

from .find_diag_var import find_diag_var, GroupedVariant, _parse_group_data
from .print_align import render_variant

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
    def sliding_window(cls, variants, groups, reference, span, flank=100):
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

        def increment(region):
            # Move one variant to spacer and add next to upstream queue
            region.variants.append(region.upstream.popleft())
            # Move from spacer to downstream queue once spacer is too large
            while region.region_length() > span:
                region.downstream.appendleft(region.variants.popleft())

        # Initialize sliding windows
        windows = {}
        for group in groups:
            windows[group] = cls(variants=[], group=group, reference=reference)
        # Read through variants and put into windows
        for index, variant in enumerate(variants):
            for group in groups:
                windows[group].upstream.append(variant)
                if index + 1 >= flank:  # Once the upstream is full
                    increment(windows[group])
                    if len(windows[group].variants) > 0:  # Can be 0 when a single var is larger than the span
                        yield cls(variants=windows[group].variants, group=group, reference=reference,
                                  upstream=windows[group].upstream, downstream=windows[group].downstream)
        # Clear the upstream window
        for index in range(flank - 1):
            for group in groups:
                increment(windows[group])
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
        var_diffs = [l - v.max_allele_len(group) for v, l in zip(vars_in_range, var_lens)]
        vars_in_range = [x for _, x in sorted(zip(var_ends, vars_in_range), key=lambda pair: pair[0])]

        # Adjust start and end positions of sequence to span all variants
        seq_ref_start = min(var_starts + [start])
        seq_ref_end = max(var_ends + [end])

        # Apply variants to reference sequence in reverse order
        out_seq = list(reference[chrom].seq[seq_ref_start:seq_ref_end+1].lower())
        for var in reversed(vars_in_range):
            replace_start = var.variant.pos - 1 - seq_ref_start
            replace_end = replace_start + len(var.variant.ref)
            if group is None:
                consensus = var.variant.ref
            else:
                alleles = var.allele_counts[group]
                if len(alleles) == 0:  # If no data for this position, use N of equal length to reference
                    consensus = "N" * var.variant.rlen
                else:
                    consensus = collapse_to_iupac(alleles.keys())
                if annotate:
                    if var.diagnostic[group] is None:
                        replacement = consensus.upper()
                    else:
                        replacement = '<' + ";".join([k+str(v) for k, v in alleles.items()]) + '>'
                else:
                    replacement = consensus.upper()
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


def run_primer3(template, target_start, target_len, options=None):
    if options is None:
        global_options = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1, # 1 == True
            # 'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_LIBERAL_BASE': 1,
            'PRIMER_OPT_SIZE': 30,
            'PRIMER_MIN_SIZE': 25,
            'PRIMER_MAX_SIZE': 35,
            # 'PRIMER_INTERNAL_MAX_SIZE': len(crrna_seq),
            'PRIMER_OPT_TM': 60.5,
            'PRIMER_MIN_TM': 53.0,
            'PRIMER_MAX_TM': 68.0,
            'PRIMER_MIN_GC': 30.0,
            'PRIMER_MAX_GC': 70.0,  # less than 50% idealy
            'PRIMER_MAX_POLY_X': 3,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
            'PRIMER_MAX_SELF_ANY_TH': 40,  # describes the tendency of a primer to bind to itself
            'PRIMER_MAX_SELF_END_TH': 40,  # 3' binding to itself (primer dimers)
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': 40,  # the tendency of the left primer to bind to the right primer
            'PRIMER_PAIR_MAX_COMPL_END_TH': 40,  # primer heterodimers
            'PRIMER_MAX_HAIRPIN_TH': 40,  # hairpins
            'PRIMER_PRODUCT_SIZE_RANGE': [[70, 150]],
        }
    else:
        global_options = parse_primer3_settings(options)

    p3_output = primer3.bindings.designPrimers(
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
    ref_diff_offset = 0  # cumulative length differences between reference and variant alleles
    var_pos_diff = 0
    for nearby_var in nearby_vars:
        # Check if variant is further away than the maximum offset
        var_pos_diff = abs(border_var.variant.pos - nearby_var.variant.pos)
        if var_pos_diff + ref_diff_offset >= max_offset:
            return {"ref": max_offset - ref_diff_offset, "group": max_offset}
        # If variant is not conserved, return pos right before variant
        if nearby_var.conserved[group] is None:
            return {"ref": var_pos_diff - len(nearby_var.variant.ref) - ref_diff_offset,
                    "group": var_pos_diff - len(nearby_var.variant.ref)}
        # Keep track of differences in allele length between the group and the reference
        max_allele_len = nearby_var.max_allele_len(group)
        ref_diff_offset += max_allele_len - len(nearby_var.variant.ref)
    return {"ref": var_pos_diff - ref_diff_offset, "group": var_pos_diff}


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

    def __init__(self,  variants, group, reference, upstream, downstream, p3, crrna_seq, downstream_seq, upstream_seq, temp_range, crrna_range):
        super().__init__(variants, group, reference, upstream, downstream)
        self.p3 = p3
        self.downstream_seq = downstream_seq
        self.crrna_seq = crrna_seq
        self.upstream_seq = upstream_seq
        self.temp_range = temp_range
        self.crrna_range = crrna_range
        self.type = "Diagnostic"

    @staticmethod
    def from_grouped_region(region, p3, crrna_seq, downstream_seq, upstream_seq, temp_range, crrna_range):
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
                               crrna_range=crrna_range)

    def left_range(self):
        """Start/stop reference sequence indexes of the left primer."""
        start = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=self.p3['PRIMER_LEFT_0'][0])
        end = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0], offset=sum(self.p3['PRIMER_LEFT_0']) - 1)
        return [start, end]

    def right_range(self):
        """Start/stop reference sequence indexes of the right primer."""
        start = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0],
                                               offset=self.p3['PRIMER_RIGHT_0'][0] - self.p3['PRIMER_RIGHT_0'][1])
        end = self.ref_pos_from_group_offset(ref_pos=self.temp_range[0],
                                             offset=self.p3['PRIMER_RIGHT_0'][0] + 1)
        return [start, end]

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
                     min_groups=1,
                     min_samples=3,
                     min_reads=5,
                     min_geno_qual=30,
                     min_map_qual=50,
                     min_freq=0.95,
                     spacer_len=28,
                     min_amp_size=80,
                     max_amp_size=300,
                     snp_offset=2,
                     offset_left=0,
                     offset_right=2):
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
    min_groups : int, optional
        The minimum number of `groups` that must be uniquely distinguished by
        the same cluster of variants.
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
    min_freq : float, optional
        The proportion of reads a variant must have to be considered real.
        This is meant to counter sequencing errors.
    spacer_len: int, optional
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
    window_width = spacer_len - offset_right - offset_left
    vcf_reader = GroupedVariant.from_vcf(variants, groups=groups,
                                         min_samples=min_samples,
                                         min_reads=min_reads,
                                         min_geno_qual=min_geno_qual)
    windower = GroupedRegion.sliding_window(vcf_reader,
                                            groups=groups.keys(),
                                            reference=reference,
                                            span=window_width,
                                            flank=max_amp_size)
    for region in windower:

        # Are there enough diagnostic variants?
        n_diag_var = sum([x is not None for x in region.diagnostic()])
        if n_diag_var < min_vars:
            region.type = 'Undiagnostic'
            yield region
            continue

        # Are all the variants in the spacer conserved?
        if any([x is None for x in region.conserved()]):
            region.type = 'Unconserved crRNA'
            yield region
            continue

        # # Is there conserved sequence adjacent to the variants to fit the crRNA?
        # region_len = region.region_length()
        # overhang_left = spacer_len - offset_right - region_len
        # if not is_nearby_conserved(region.group, region.variants[-1], region.upstream, offset_right + 1):
        #      region.type = 'Unconserved crRNA'
        #      yield region
        # if not is_nearby_conserved(region.group, region.variants[0], region.downstream, overhang_left + 1):
        #      region.type = 'Unconserved crRNA'
        #      yield region

        # Is there conserved adjacent sequence for the crRNA?
        overhang_left = spacer_len - region.region_length() - offset_right
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
            region.type = 'Unconserved crRNA'
            yield region
            continue

        # Is there enough adjacent conserved regions to design primers?
        consv_len_up = consv_border_n(group=region.group,
                                               border_var=region.variants[-1],
                                               nearby_vars=region.upstream,
                                               max_offset=max_amp_size)
        consv_len_dn = consv_border_n(group=region.group,
                                               border_var=region.variants[0],
                                               nearby_vars=region.downstream,
                                               max_offset=max_amp_size)
        if consv_len_up["group"] - overhang_len_up['group'] < 30:  #TODO base on primer3 primer size
            region.type = 'Unconserved seq'
            yield region
            continue
        if consv_len_dn["group"] - overhang_len_dn['group'] < 30:  #TODO base on primer3 primer size parameters
            region.type = 'Unconserved seq'
            yield region
            continue

        # Run primer3 on group-specific template
        start_crrna_ref = region.variants[0].variant.pos - 1 - overhang_len_dn['ref']
        end_crrna_ref = region.variants[-1].variant.pos - 1 + overhang_len_up['ref']
        start_tmp_ref = region.variants[0].variant.pos - 1 - consv_len_dn['ref']
        end_tmp_ref = region.variants[-1].variant.pos - 1 + consv_len_up['ref']
        downstream_seq = region.sequence(reference=reference, start=start_tmp_ref, end=start_crrna_ref-1, group=region.group)
        upstream_seq = region.sequence(reference=reference, start=end_crrna_ref+1, end=end_tmp_ref, group=region.group)
        crrna_seq = region.sequence(reference=reference, start=start_crrna_ref, end=end_crrna_ref, group=region.group)
        template_seq = downstream_seq + crrna_seq + upstream_seq
        start_crrna_tmp = len(downstream_seq)
        end_crrna_tmp = start_crrna_tmp + len(crrna_seq) - 1

        p3_out = run_primer3(template_seq, target_start=start_crrna_tmp, target_len=len(crrna_seq))
        if p3_out['PRIMER_PAIR_NUM_RETURNED'] == 0:
            region.type = 'No primers'
            yield region
            continue

        # Return data for the diagnostic region
        region.type = 'Diagnostic'
        output = DiagosticRegion.from_grouped_region(region, p3=p3_out, crrna_seq=crrna_seq, downstream_seq=downstream_seq, upstream_seq=upstream_seq,
                                                     temp_range=[start_tmp_ref, end_tmp_ref], crrna_range=[start_crrna_ref, end_crrna_ref])

        yield output


def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Find regions where there are conserved variants for each group that are not found in other groups.')
    parser.add_argument('metadata', type=str,
                        help='A TSV file containing data with one row per sample. Two columns are required: `sample_id`, which contains the same sample IDs used in the VCF file, and `group`, which identifies which group each sample belongs to.')
    parser.add_argument('vcfs', type=str, nargs="+",
                        help='One or more VCF files containing variant data for the samples grouped by the `metadata` file.')
    parser.add_argument('--groups', type=str, nargs="+",
                        help='One or more groups that are to be distinguished by variants. These should match the values of the `group` column in the `metadata` file.')
    parser.add_argument('--reference', type=str,
                        help='The reference file used to make the `vcfs` VCF file. If supplied, the sequence of the region containing each conserved variant is returned with the output.')
    parser.add_argument('--out', type=str,
                        help='The output file to create. (default: print to screen)')
    parser.add_argument('--log', type=str,
                        help='The location to save the log file (the contents of the standard error stream). (default: print to screen)')
    return parser.parse_args()


def parse_vcfs(paths):
    """Create a dict of file handles for a list of paths to VCF files"""
    out = {}
    for path in paths:
        out[path] = pysam.VariantFile(path)
    return out


def _format_p3_output(p3_out):
    """Reformat data for best primer pair for CSV output"""
    return {primer3_col_key[n]: p3_out[n] for n in primer3_col_names}


def _render_header(stream, out_sep=','):
    print("group", "chrom", "n_diag", "fwd_start", "fwd_end", "rev_start", "rev_end", "crrna_start", "crrna_end", "seq_start", "seq_end",
          *[primer3_col_key[n] for n in primer3_col_names],
          sep=out_sep, file=stream)


@contextmanager
def writer(file_path=None, default_stream=sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "w")
    yield file_handle
    if file_path is not None:
        file_handle.close()


def _format_for_csv(region, reference):
    fwd_range = region.left_range()
    rev_range = region.right_range()
    crrna_range = region.crrna_range
    temp_range = region.temp_range

    group = region.group
    chrom = region.variants[0].variant.chrom
    n_diag = sum([x is not None for x in region.diagnostic()])

    fwd_start = fwd_range[0] + 1  # The +1 goes from 0-based in 1-based indexing
    fwd_end = fwd_range[1] + 1
    rev_start = rev_range[0] + 1
    rev_end = rev_range[1] + 1
    crrna_start = crrna_range[0] + 1
    crrna_end = crrna_range[1] + 1
    seq_start = region.temp_range[0] + 1
    seq_end = region.temp_range[1] + 1

    # up_seq = ''.join(region.upstream_seq)
    # down_seq = ''.join(region.downstream_seq)
    # crrna_seq = ''.join(region.crrna_seq)

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

    # Modify these values to change columns and their names:
    output = {
        "group": group,
        "chrom": chrom,
        "n_diag": n_diag,
        "reg_from": seq_start,
        "reg_to": seq_end,
        "diag_from": crrna_start,
        "diag_to": crrna_end,
        "fwd_from": fwd_start,
        "fwd_to": fwd_end,
        "rev_from": rev_start,
        "rev_to": rev_end,
        "seq_adj_left": seq_adj_left,
        "seq_primer_fwd": seq_primer_left,
        "seq_inter_left": seq_amp_left,
        "seq_diag": seq_crrna,
        "seq_inter_right": seq_amp_right,
        "seq_primer_rev": seq_primer_right,
        "seq_adj_right": seq_adj_right
    }
    output.update(_format_p3_output(region.p3))

    # if seq_start == 42133:
    #     import pdb; pdb.set_trace()


    return output


def run_all():
    # Read command line arguments
    args = parse_command_line_args()

    # Prepare input data
    vcf_handles = parse_vcfs(args.vcfs)
    first_vcf_handle = pysam.VariantFile(args.vcfs[0])
    first_var = next(first_vcf_handle)
    groups = _parse_group_data(args.metadata, groups=args.groups, possible=first_var.samples.keys())
    reference = _parse_reference(args.reference)
    stats = defaultdict(int)

    # Prepare output
    with writer(args.out, sys.stdout) as output_stream, writer(args.log, sys.stderr) as log_stream:
        for vcf_handle in vcf_handles.values():
            for region in find_diag_region(vcf_handle, groups, reference=reference):
                stats[region.type] += 1
                print(stats)
                if region.type == "Diagnostic":
                    # seqs = {g: region.sequence(reference=reference, start=region.temp_range[0], end=region.temp_range[1], group=g) for g in groups.keys()}
                    # ref = region.sequence(reference=reference, start=region.temp_range[0], end=region.temp_range[1], group=None)
                    # render_variant(seqs, ref)

                    columns = _format_for_csv(region, reference)
                    if stats["Diagnostic"] == 1: # If first output printed
                        print(*columns.keys(), sep=',', file=output_stream)
                    print(*columns.values(), sep=',', file=output_stream)


def main():
    # import cProfile; cProfile.runctx('run_all()', globals=globals(), locals=locals(), sort='cumulative')
    run_all()


if __name__ == "__main__":
    # Test command:
    # python -m diagvar.find_diag_region test_data/test_metadata.tsv test_data/unfilt_allscafs_n666.vcf.gz --groups NA1 NA2 EU1 EU2 --reference test_data/PR-102_v3.1_s0001.fasta --out test.csv
    main()

# if __name__ == "__main__":
#     vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
#     groups = _read_group_data('test_data/test_metadata.tsv')
#     groups = {g: v for g, v in groups.items() if g in ["NA1", "NA2", "EU1", "EU2"]}
#     ref = 'test_data/PR-102_v3.1.fasta'
#     find_diag_region(vcf,
#                      groups=groups,
#                      reference=ref)

