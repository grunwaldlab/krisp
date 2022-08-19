#! /bin/env python3

import sys
import argparse
import pandas
import pysam
from contextlib import contextmanager
from collections import Counter
from .find_diag_region import find_diag_region

# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM = ('(', ')')
CRRNA_DELIM = ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'


def _parse_command_line_args():
    parser = argparse.ArgumentParser(
        description='Find conserved variants for each group that are not '
                    'found in other groups.'
    )
    parser.add_argument(
        'metadata', type=str,
        help='A TSV file containing data with one row per sample. '
             'Two columns are required: `sample_id`, which contains the '
             'same sample IDs used in the VCF file, and `group`, which '
             'identifies which group each sample belongs to.'
    )
    parser.add_argument(
        'vcf', type=str, nargs="+",
        help='One or more VCF files containing variant data for the samples '
             'grouped by the `metadata` file.'
    )
    parser.add_argument(
        '--groups', type=str, nargs="+",
        help='One or more groups that are to be distinguished by variants. '
             'These should match the values of the `group` column in the '
             '`metadata` file. (default: use all groups in the `metadata` '
             'file)'
    )
    parser.add_argument(
        '--reference', type=str, 
        help='The reference file associated with the `VCF` file. '
             'If supplied, the sequence of the region containing each '
             'conserved variant is returned with the output.'
    )
    parser.add_argument(
        '--out', type=str,
        help='The output file to print to. (default: print to screen)'
    )
    parser.add_argument(
        '--log', type=str,
        help='The file to print error messages and status updates to '
             '(default: print to screen via standard error)'
    )
    parser.add_argument(
        '--primer3', action=argparse.BooleanOptionalAction, default=False,
        help='Use Primer3 to pick primers to amplify each variant. '
             'Only variants for which Primer3 can find primers will '
             'be included in the output.'
    )
    parser.add_argument(
        '--primer3_options', type=str,
        help='The path to a file with settings for Primer3.'
    )
    parser.add_argument(
        '--flank', type=int, default=100,
        help='The length of the region on either side of the variant '
             'that must be conserved within the group. (default: %(default)s)'
    )
    parser.add_argument(
        '--min_vars', type=int, default=1,
        help='The number of diagnostic variants that must appear in the '
             'spacer sequence. (default: %(default)s)'
    )
    parser.add_argument(
        '--min_groups', type=int, default=1,
        help='The number of groups that must be uniquely distinguished '
             'by each variant. (default: %(default)s)'
    )
    parser.add_argument(
        '--min_samples', type=int, default=5,
        help='The number of samples with acceptable data (see `--min_reads`)'
             ' each group must have for a given variant. '
             '(default: %(default)s)'
    )
    parser.add_argument(
        '--min_reads', type=int, default=10,
        help='The number of reads a variant must be represented by in a '
             'given sample for that sample to count towards `--min_samples`. '
             'This corresponds to the per-sample GATK output in the VCF '
             'encoded as "DP". (default: %(default)s)'
    )
    parser.add_argument(
        '--min_geno_qual', type=int, default=40,
        help='The minimum genotype quality score (phred scale). This '
             'corresponds to the per-sample GATK output in the VCF encoded as '
             '"GQ". (default: %(default)s)'
    )
    parser.add_argument(
        '--min_map_qual', type=int, default=50,
        help='The minimum root square mean mapping quality of reads '
             'supporting a variant (phred scale). This corresponds to '
             'the per-variant GATK output in the VCF encoded as "MQ". '
             '(default: %(default)s)'
    )
    parser.add_argument(
        '--min_freq', type=float, default=0.05,
        help='The proportion of reads a variant must have to be considered '
             'real in a given sample. (default: %(default)s)'
    )
    parser.add_argument(
        '--spacer_len', type=int, default=28,
        help='The length of the spacer sequence that the crRNA will bind to. '
             '(default: %(default)s)'
    )
    parser.add_argument(
        '--snp_offset', type=int, default=-3,
        help='How many bases from the 5` end of the spacer will the '
             'SNP occur. Negative values will index from the 3` end. '
             'For example 0 is the 5` end of the spacer and -1 is the 3`  '
             'end. (default: %(default)s)'
    )
    parser.add_argument(
        '--annotation', default=['default'], nargs=1,
        choices=['none', 'minimal', 'default', 'all'],
        help='What type of annotations are used in the sequence output. '
             '"none": The sequence is returned is the best guess consensus '
             'sequence for the target group. '
             '"minimal": oligos are marked and diagnostic variant is '
             'uppercase. '
             '"default": oligos are marked, unambiguous variants are '
             'uppercase, ambiguous variants are annotated with the '
             'counts for each allele. '
             '"all": oligos are marked, all variants are annotated with '
             'the counts for each allele.'
    )
    return parser.parse_args()


def _validate_command_line_args(original_args):
    return original_args


#   https://stackoverflow.com/questions/22264504/using-python-with-statement-with-sys-stdout
@contextmanager
def writer(file_path=None, default_stream=sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "w")
    yield file_handle
    if file_path is not None:
        file_handle.close()


def _print_stats_header(stream, args, variant_stat_names):
    max_nchar = max([len(n) for n in variant_stat_names + args.groups])
    header_parts = [n.ljust(max_nchar) for n in variant_stat_names +
                    args.groups]
    print('| '.join(header_parts), file=stream)


def _print_stdout_header(args, stream, primer3_col_key, sep='\t'):
    col_names = ["chromosome", "position", "target", *args.groups]
    if args.reference is not None:
        col_names += [*args.groups, "reference", "spacer"]
    if args.primer3:
        col_names += [primer3_col_key[n] for n in primer3_col_key.values()]
    print(col_names, sep=sep, file=stream)


def _print_log_stats(stream, variant_counts, group_counts, endline=False):
    max_nchar = max([len(n) for n in variant_counts.keys() +
                     group_counts.keys()])
    var_info = [str(x).ljust(max_nchar) for x in variant_counts.values()]
    group_info = [str(x).ljust(max_nchar) for x in group_counts.values()]
    print('| '.join(var_info + group_info), file=stream,
          end='\n' if endline else '\r')
  
  
def _get_group_map(metadata_path, sample_col="sample_id", group_col="group"):
    """
    Reads metadata file and returns dictionary with sample IDs as keys and group
    IDs as items
    """
    metadata = pandas.read_csv(metadata_path, sep='\t')
    return dict(zip(metadata[sample_col].tolist(),
                    metadata[group_col].tolist()))


def _print_result(result):
    pass


def _command_line_entry_point():
    """What is called when this module is executed as a script."""
    # Parse and verify the command line input.
    original_args = _parse_command_line_args()
    args = _validate_command_line_args(original_args)
    breakpoint()
    # Initialize counters for statistics printed to stderr or the log file.
    variant_stat_names = [
        'Processed',
        'Undiagnostic',
        'Missing data',
        'Low map qual',
        'Unconserved',
        'No primers',
    ]
    variant_counts = Counter({k: 0 for k in variant_stat_names})
    group_map = _get_group_map(args.metadata)
    group_counts = Counter(group_map.values())
    
    # Set up primer3 columns to print
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
        'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH'
    ]
    primer3_col_key = {n: n.replace("PRIMER_", "").replace("_0", "").lower()
                       for n in primer3_col_names}
  
    with writer(args.out, sys.stdout) as output_stream, \
            writer(args.log, sys.stderr) as log_stream:
        # print header if needed
        _print_stdout_header(output_stream, args, primer3_col_key)
        if args.out is not None and args.log is None:
            _print_stats_header(log_stream, args, variant_stat_names)
        # Process each input file
        for vcf_file in args.vcfs:
            vcf_handle = pysam.VariantFile(vcf_file)
            for result in find_diag_region(vcf_handle, group_map,
                                           reference=args.reference,
                                           primer3=args.primer3,
                                           min_vars=args.min_vars,
                                           min_groups=args.min_groups,
                                           min_samples=args.min_samples,
                                           min_reads=args.min_reads,
                                           min_geno_qual=args.min_geno_qual,
                                           min_map_qual=args.min_map_qual,
                                           min_freq=args.min_freq,
                                           spacer_len=args.spacer_len,
                                           snp_offset=args.snp_offset):
                _print_result(result)
        if args.out is None or args.log is not None:
            _print_stats_header(log_stream, args, variant_stat_names)
        _print_log_stats(log_stream, variant_counts, group_counts, endline=True)


if __name__ == "__main__":
    _command_line_entry_point()
