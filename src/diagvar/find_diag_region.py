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
PRIMER_DELIM = ('(', ')')
CRRNA_DELIM = ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'
iupac_key = {tuple((x for x in sorted(v))): k for k, v in
             IUPACData.ambiguous_dna_values.items()}
iupac_key[(UNKNOWN_CHAR,)] = 'N'


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

    @classmethod
    def sliding_window(cls, variants, groups, reference, span, flank=100):
        """Generate GroupedRegion objects as a sliding window along variants

        Parameters
        ----------
        variants : iterator returning GroupedVariant
            The variants to generate GroupedRegion objects from.
        groups : list of str
            The groups to return sliding windows for
        reference : dict
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
                    yield windows[group]
        # Clear the upstream window
        for index in range(flank - 1):
            for group in groups:
                increment(windows[group])
                yield windows[group]

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
        # Get length of the reference spanned by these variants
        var_poses = [x.variant.pos for x in self.variants]
        out = max(var_poses) - min(var_poses) + 1

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
                ref_seq = ref_seq[:replace_start] + [rend] + \
                          ref_seq[replace_end:]
            return ref_seq

    def sequence(self, reference, start, end, group=None, var_upper=True):
        """Infer the sequence for the current group based on the reference sequence and variants
        """
        if group is None:
            group = self.group
        # Find variants within target region
        all_vars = self.downstream + self.variants + self.upstream
        vars_in_range = [x for x in all_vars if start <= x.variant.pos <= end]
        # Check that variants are all on the same chormosome
        if len({x.variant.chrom for x in vars_in_range}) > 1:
            raise ValueError('Variants cannot span multiple chromosomes')
        chrom = self.variants[-1].variant.chrom
        # Insert consensus for each variant in reference sequence
        ref_seq = reference[chrom].seq[start:end].lower()
        if var_upper:
            ref_seq = list(ref_seq.lower())
        else:
            ref_seq = list(ref_seq)
        for var in vars_in_range:
            replace_start = var.variant.pos - start
            replace_end = replace_start + len(var.variant.ref)
            consensus = collapse_to_iupac(var.allele_counts[group].keys())
            if var_upper:
                consensus = consensus.upper()
            ref_seq = ref_seq[:replace_start] + [consensus] + ref_seq[replace_end:]
        return ref_seq

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


def _check_variant_cluster(variants, subset):
    """Performs a series of checks to see if the cluster is diagnostic"""

    pass


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
        options = dict([l.strip().split('=') for l in handle.readlines()])
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


def run_primer3(template, target_start, target_end, options=None):
    if options is None:
        global_options = {
            'PRIMER_LIBERAL_BASE': 1,
            'PRIMER_OPT_SIZE': 30,
            'PRIMER_MIN_SIZE': 25,
            'PRIMER_MAX_SIZE': 35,
            'PRIMER_OPT_TM': 60.5,
            'PRIMER_MIN_TM': 54.0,
            'PRIMER_MAX_TM': 67.0,
            'PRIMER_MIN_GC': 30.0,
            'PRIMER_MAX_GC': 50.0,  # less than 50% idealy
            'PRIMER_MAX_POLY_X': 3,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
            'PRIMER_MAX_SELF_ANY_TH': 30,  # describes the tendency of a primer to bind to itself
            'PRIMER_MAX_SELF_END_TH': 30,  # 3' binding to itself (primer dimers)
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': 30,  # the tendency of the left primer to bind to the right primer
            'PRIMER_PAIR_MAX_COMPL_END_TH': 30,  # primer heterodimers
            'PRIMER_MAX_HAIRPIN_TH': 30,  # hairpins
            'PRIMER_PRODUCT_SIZE_RANGE': [[80, 140]]
        }
    else:
        global_options = parse_primer3_settings(options)

    p3_output = primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': "".join(template),
            'SEQUENCE_TARGET': [target_start, target_end - target_start + 1]
        },
        global_options
    )
    return p3_output


def consv_ref_border_pos(group, border_var, nearby_vars, max_offset):
    """Find the maximum size of adjacent conserved reference sequence"""
    ref_diff_offset = 0  # cumulative differences between reference and variant alleles
    for nearby_var in nearby_vars:
        # Check if variant is further away than the maximum offset
        var_pos_diff = abs(border_var.variant.pos - nearby_var.variant.pos)
        if var_pos_diff + ref_diff_offset >= max_offset:
            return max_offset - ref_diff_offset
        # If variant is not conserved, return pos right before variant
        if nearby_var.conserved[group] is None:
            return var_pos_diff - len(nearby_var.variant.ref)
        # Keep track of differences in allele length between the group and the reference
        max_allele_len = nearby_var.max_allele_len(group)
        ref_diff_offset += max_allele_len - len(nearby_var.variant.ref)


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
                     snp_offset=2,
                     offset_left=0,
                     offset_right=3):
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
    reference : str
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
    flank = 100  # TODO: base on max amplicon size
    max_amplicon_len = 1000  # TODO: base on max amplicon size
    reference = _parse_reference(reference)
    window_width = spacer_len - offset_right - offset_left
    vcf_reader = GroupedVariant.from_vcf(variants, groups=groups,
                                         min_samples=min_samples,
                                         min_reads=min_reads,
                                         min_geno_qual=min_geno_qual)
    windower = GroupedRegion.sliding_window(vcf_reader,
                                            groups=groups.keys(),
                                            reference=reference,
                                            span=window_width,
                                            flank=flank)
    for region in windower:

        # Are there enough diagnostic variants?
        n_diag_var = sum([x is not None for x in region.diagnostic()])
        if n_diag_var < min_vars:
            continue

        # Are all the variants in the spacer conserved?
        if any([x is None for x in region.conserved()]):
            continue

        # Is there conserved sequence adjacent to the variants to fit the crRNA?
        region_len = region.region_length()
        overhang_left = spacer_len - offset_right - region_len
        if not is_nearby_conserved(region.group, region.variants[-1], region.upstream, offset_right + 1):
            continue
        if not is_nearby_conserved(region.group, region.variants[0], region.downstream, overhang_left + 1):
            continue

        # Can primers be designed in adjacent conserved regions?
        num_bp_consv_up = consv_ref_border_pos(group=region.group,
                                               border_var=region.variants[-1],
                                               nearby_vars=region.upstream,
                                               max_offset=max_amplicon_len)
        num_bp_consv_dn = consv_ref_border_pos(group=region.group,
                                               border_var=region.variants[0],
                                               nearby_vars=region.downstream,
                                               max_offset=max_amplicon_len)
        spacer_len_diff = sum([x.max_allele_len(region.group) - len(x.variant.ref) for x in region.variants])
        spacer_var_len = region.region_length()
        overhang_left = spacer_len - spacer_var_len - offset_right
        num_bp_overhang_up = consv_ref_border_pos(border_var=region.variants[-1],
                                               nearby_vars=region.upstream,
                                               max_offset=offset_right)
        num_bp_overhang_dn = consv_ref_border_pos(border_var=region.variants[0],
                                               nearby_vars=region.downstream,
                                               max_offset=overhang_left)

        # Get location of crrna and template sequence in terms of reference
        start_crrna_ref = region.variants[0].variant.pos - num_bp_overhang_dn
        end_crrna_ref = region.variants[-1].variant.pos + num_bp_overhang_up
        start_tmp_ref = region.variants[0].variant.pos - num_bp_consv_dn
        end_tmp_ref = region.variants[-1].variant.pos + num_bp_consv_up
        print(start_crrna_ref, end_crrna_ref, start_tmp_ref, end_tmp_ref)

        # Get template sequence for downstream, crrna, and upstream
        upstream_seq = region.sequence(reference=reference, start=start_tmp_ref, end=start_crrna_ref-1)
        downstream_seq = region.sequence(reference=reference, start=end_crrna_ref+1, end=end_tmp_ref)
        crrna_seq = region.sequence(reference=reference, start=start_crrna_ref, end=end_crrna_ref)
        template_seq = downstream_seq + crrna_seq + upstream_seq

        # Run primer3 on group-specific template
        start_crrna_tmp = len("".join(downstream_seq))
        end_crrna_tmp = start_crrna_tmp + len("".join(crrna_seq)) - 1
        p3_out = run_primer3(template_seq, start_crrna_tmp, end_crrna_tmp)
        if p3_out['PRIMER_PAIR_NUM_RETURNED'] == 0:
            continue


        print(crrna_seq)
