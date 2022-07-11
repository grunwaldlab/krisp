# External dependencies
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

# Command line interface
parser = argparse.ArgumentParser(description='Find conserved variants for each group that are not found in other groups.')
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
parser.add_argument('--primer3', action=argparse.BooleanOptionalAction, default=False,
                    help='Use Primer3 to pick primers to amplify each variant. Only variants for which Primer3 can find primers will be included in the output.')
parser.add_argument('--primer3_options', type=str,
                    help='The path to a file with settings for Primer3.')
parser.add_argument('--format',
                    default='tsv',
                    nargs=1,
                    choices=['tsv', 'csv', 'vcf', 'human'],
                    help='Which format to print the output. vcf: filter input VCF to only SNPs that distinguish groups. human: an information rich format meant for maunal interpretation.')
parser.add_argument('--flank', type=int, default=100,
                    help='The length of the region that must be conserved within the group, on either side of the group-specific variant. (default: %(default)s)')
parser.add_argument('--min_groups', type=int, default=1,
                    help='The number of groups that must uniquly distinguished by a location for the variant to be recorded. (default: %(default)s)')
parser.add_argument('--min_samples', type=int, default=5,
                    help='The number of samples with acceptable data (see `--min_reads`) each group must have for a given variant. (default: %(default)s)')
parser.add_argument('--min_reads', type=int, default=10,
                    help='The number of reads a variant must be represented by in a given sample for the data of that sample to be considered. This corresponds to the per-sample GATK output in the VCF encoded as "DP". (default: %(default)s)')
parser.add_argument('--min_geno_qual', type=int, default=40,
                    help='The minimum genotype quality score (phred scale). This corresponds to the per-sample GATK output in the VCF encoded as "GQ". (default: %(default)s)')
parser.add_argument('--min_map_qual', type=int, default=50,
                    help='The minimum root square mean mapping quality of reads supporting a variant (phred scale). This corresponds to the per-variant GATK output in the VCF encoded as "MQ". (default: %(default)s)')
parser.add_argument('--min_freq', type=float, default=0.05,
                    help='The proportion of reads a variant must have to be considered real. This is meant to counter sequencing errors. (default: %(default)s)')
parser.add_argument('--spacer_len', type=int, default=28,
                    help='The length of the spacer sequence that the crRNA will bind to and the primers will amplify. (default: %(default)s)')
parser.add_argument('--snp_offset', type=int, default=-3,
                    help='How many bases from the 5` end of the spacer will the SNP occur. For example 0 is the 5` end of the spacer and -1 is the 3` end. (default: %(default)s)')
parser.add_argument('--annotation',
                    default=['default'],
                    nargs=1,
                    choices=['none', 'minimal', 'default', 'all'],
                    help='What type of annotations are used in the sequence output. "none": The sequence is returned is the best guess consensus sequence for the target group. "minimal": oligos are marked and diagnositic variant is uppercase. "default": oligos are marked, unambigous variants are uppercase, ambiguous variants are annotated with the counts for each allele. "all": oligos are marked, all variants are annotated with the counts for each allele.')
          
          
# Constants
SNP_DELIM = ('<', '>')
PRIMER_DELIM =  ('(', ')')
CRRNA_DELIM =  ('{', '}')
HETERO_DELIM = '/'
UNKNOWN_CHAR = '?'

def to_number_if_can(x):
  try:
    if int(float(x)) == float(x) and '.' not in x:
      return int(x)
    else:
      return float(x)
  except ValueError:
    return x

def parse_primer3_settings(file):
  """
  Reads primer3 BoulderIO format, assuming only global settings (starts with PRIMER_),
  and returns a dict in the format used by primer3-py.
  """
  with open(args.primer3_options) as handle:
    options = dict([l.strip().split('=') for l in handle.readlines()])
    for opt, val in options.items():
      if ' ' in val or  ';' in val:
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

def run_primer3(template, site_index, site_len):
  if args.primer3_options is None:
    global_options = {
      'PRIMER_LIBERAL_BASE': 1,
      'PRIMER_OPT_SIZE': 30,
      'PRIMER_MIN_SIZE': 25,
      'PRIMER_MAX_SIZE': 35,
      'PRIMER_OPT_TM': 60.5,
      'PRIMER_MIN_TM': 54.0,
      'PRIMER_MAX_TM': 67.0,
      'PRIMER_MIN_GC': 30.0,
      'PRIMER_MAX_GC': 50.0, # less than 50% idealy
      'PRIMER_MAX_POLY_X': 3,
      'PRIMER_MAX_NS_ACCEPTED': 0,
      'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
      'PRIMER_MAX_SELF_ANY_TH': 30, # describes the tendency of a primer to bind to itself
      'PRIMER_MAX_SELF_END_TH': 30, # 3' binding to itself (primer dimers)
      'PRIMER_PAIR_MAX_COMPL_ANY_TH': 30, # the tendency of the left primer to bind to the right primer
      'PRIMER_PAIR_MAX_COMPL_END_TH': 30, # primer heterodimers
      'PRIMER_MAX_HAIRPIN_TH': 30,  # hairpins
      'PRIMER_PRODUCT_SIZE_RANGE': [[80,140]]
    }
  else:
    global_options = parse_primer3_settings(args.primer3_options)
  
  p3_output = primer3.bindings.designPrimers(
    {
      'SEQUENCE_TEMPLATE': str(template),
      'SEQUENCE_TARGET': [site_index, site_len]
    },
    global_options
    )
  return p3_output


def get_group_map(metadata_path, sample_col="sample_id", group_col="group"):
  """
  Reads metadata file and returns dictionary with sample IDs as keys and group IDs as items
  """
  metadata = pandas.read_csv(metadata_path, sep='\t')
  return dict(zip(metadata[sample_col].tolist(), metadata[group_col].tolist()))

def get_reference(ref_path):
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



def group_genotypes(record, group_map):
  """
  For a given variant, return a dictionary with the counts of each genotype for each group.
  """
  genotypes = {g: {} for g in set(group_map.values())}
  for sample_id, data in record.samples.items():
    if data['DP'] == 0: # https://gatk.broadinstitute.org/hc/en-us/articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot
      alleles = UNKNOWN_CHAR 
    else:
      alleles = sorted(list(set(data.alleles)))
      alleles = [UNKNOWN_CHAR if a is None else a for a in alleles]
      alleles = "/".join(alleles) # heterzygotes will have a / in their genotype
    if alleles in genotypes[group_map[sample_id]]:
      genotypes[group_map[sample_id]][alleles] += 1
    else:
      genotypes[group_map[sample_id]][alleles] = 1
  return genotypes



def extract_seq(record, group, downstream, upstream, annotate = 'default'):
  """
  Return the inferred consensus sequence for the group based on the reference sequence and surrouding SNPs.
  Lower case letters correpond to bases matching the reference and upper case letters are differences.
  The output is a list with one element for each position in the reference
  """
  # Get reference sequence surrounding variant
  sequence = list(reference[record.chrom].seq[record.pos - args.flank - 1 : record.pos + args.flank])
  sequence = [s.lower() for s in sequence]
  
  # Modify the reference sequences so that it matches the group consensus sequence
  iupac_key = {v:k for k,v in IUPACData.ambiguous_dna_values.items()}
  vars_in_range = [x for x in downstream if abs(record.pos - x.pos) <= args.flank]
  vars_in_range.append(record)
  vars_in_range.extend([x for x in upstream if abs(record.pos - x.pos) <= args.flank])
  for var in vars_in_range:
    # Create representation of the genotype for the group will all the information
    var_gt = group_genotypes(var, group_map)
    group_gt = var_gt[group]
    encoded_gt = SNP_DELIM[0] + ";".join([k+str(v) for k, v in group_gt.items()]) + SNP_DELIM[1]
    # Create a siplified represenation 
    sorted_group_gt = sorted(group_gt.items(), key=lambda x: x[1], reverse=True)
    simplified_gt = sorted_group_gt[0][0] # Most common 
    if simplified_gt == "*":
      simplified_gt = ""
    elif simplified_gt == UNKNOWN_CHAR:
      simplified_gt = "N"
    elif HETERO_DELIM in simplified_gt:
      split_gt = simplified_gt.split(HETERO_DELIM)
      if len(set([len(x) for x in split_gt])) == 1: # all the same length 
        per_base_gt = [[x[i] for x in split_gt] for i in range(len(split_gt[0]))]
        simplified_gt = "".join([iupac_key["".join(sorted(x))] for x in per_base_gt])
      else:
        simplified_gt = split_gt[0]
        print(f"Warning: Multiple alleles of different length at the same site cannot be simplied to a single character. The first allele will be displayed.", file = log_stream)
    # Decide which represenation to use
    diagnostic = var == record # TODO: find way to identify nearby diagnsotic variants without redundant calculations
    unambiguous = len(group_gt) == 1 and list(group_gt.values())[0] >= args.min_samples and re.search('^[GACTgact]+$', list(group_gt.keys())[0])
    if annotate == 'none':
      replace_char = simplified_gt.lower()
    elif annotate == 'minimal':
      if diagnostic:
        replace_char = simplified_gt.upper()
      else:
        replace_char = simplified_gt.lower()
    elif annotate == 'default':
      if unambiguous and not diagnostic:
        replace_char = simplified_gt.upper()
      else:
        replace_char = encoded_gt
    elif annotate == 'all':
      replace_char = encoded_gt
    else:
      raise Exception('Invalid "annotate" value.')

    # Replace reference allele with group info
    replace_start = args.flank - (record.pos - var.pos)
    replace_end = replace_start + len(var.ref)
    sequence = sequence[:replace_start] + [replace_char] + sequence[replace_end:] 
  return sequence

#https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
def flatten(t):
    return [item for sublist in t for item in sublist]


def diagnostic_gts(genotypes):
  """
  Return a dict of conserved variants only present in a single group
  """
  unique_gts = {}
  for g in groups:
    unique_gts[g] = set(flatten([x.split(HETERO_DELIM) for x in genotypes[g].keys()]))
    if UNKNOWN_CHAR in unique_gts[g]: unique_gts[g].remove(UNKNOWN_CHAR)
  diag_gt = {k: v if len(v) == 1 else set() for k, v in copy.deepcopy(unique_gts).items()} # NOTE: ignores heterozygotes
  for group in groups:
    for other_group in groups:
      if other_group != group:
        diag_gt[group] -= unique_gts[other_group]
  return diag_gt


def process_one_variant(record, group_map, downstream, upstream):
  """
  Looks at a line in a VCF file and decides if it meets the criteria of a diagnositc SNP, as defined by the options
  """
  
  variant_stats['Processed'] += 1
  
  # Get genotypes for each group
  genotypes = group_genotypes(record, group_map)
  diag_gt = diagnostic_gts(genotypes)

  # Ignore heterozygous variants (NOTE: not actually doing anything, see above)
  diag_gt = {k: {x for x in v if HETERO_DELIM not in x} for k, v in diag_gt.items()}
  
  # Ignore deletions
  diag_gt = {k: {x for x in v if x != "*"} for k, v in diag_gt.items()}

  # Ignore insertions
  diag_gt = {k: {x for x in v if len(x) == len(record.ref)} for k, v in diag_gt.items()}
  
  # Ignore sites without enough canidate variants
  if sum([len(x) > 0 for x in diag_gt.values()]) < args.min_groups:
    variant_stats['Undiagnostic'] += 1
    return [None], None, [None]
  
  # Ignore sites with not enough reads/samples
  samples_per_group = Counter([group_map[x.name] for x in record.samples.values() if x['DP'] >= args.min_reads and x['GQ'] >= args.min_geno_qual])
  enough_samples_per_group = {g: g in samples_per_group and samples_per_group[g] >= args.min_samples for g in groups}
  if any([not x for x in enough_samples_per_group.values()]):
    variant_stats['Missing data'] += 1
    return [None], None, [None]
  
  # Ignore sites that are in low mapping quality regions
  if record.info['MQ'] < args.min_map_qual:
    variant_stats['Low map qual'] += 1
    return [None], None, [None]
  
  # Ignore sites that are not flanked by conserved sequence
  vars_in_range = [x for x in downstream if abs(record.pos - x.pos) <= args.flank]
  vars_in_range.extend([x for x in upstream if abs(record.pos - x.pos) <= args.flank])
  for var in vars_in_range:
    var_gt = group_genotypes(var, group_map)
    for group in groups:
      if len(diag_gt[group]) > 0:
        if sum([k != UNKNOWN_CHAR for k in var_gt[group].keys()]) > 1:
          variant_stats['Unconserved'] += 1
          return [None], None, [None]
  
  # Filter based on primer3 matches
  target_groups = [g for g in groups if len(diag_gt[g]) > 0]
  if args.primer3:
    primer3_out = {}
    for group in target_groups:
      sequence = extract_seq(record, group, downstream, upstream, annotate = "minimal")
      primer3_out[group] = run_primer3(template = "".join(sequence), site_index = get_spacer_start(), site_len = args.spacer_len)
      if primer3_out[group]['PRIMER_PAIR_NUM_RETURNED'] == 0:
        variant_stats['No primers'] += 1
        return [None], None, [None]
  else:
    primer3_out = [None] * len(target_groups)
  
  # Update counts
  for group in groups:
    if len(diag_gt[group]) > 0:
      group_counts[group] += 1
  
  return target_groups, genotypes, primer3_out

def get_spacer_start():
  if args.snp_offset >=0:
    return args.flank - args.snp_offset
  else:
    return args.flank - args.spacer_len - args.snp_offset

def print_log_stats_header(log_stream = sys.stderr):
  max_nchar = max([len(n) for n in variant_stat_names + groups])
  header_parts = [n.ljust(max_nchar) for n in variant_stat_names + groups]
  print('| '.join(header_parts), file = log_stream)

def print_log_stats(log_stream = sys.stderr, endline = False):
  max_nchar = max([len(n) for n in variant_stat_names + groups])
  var_info = [str(variant_stats[n]).ljust(max_nchar) for n in variant_stat_names]
  group_info = [str(group_counts[n]).ljust(max_nchar) for n in groups]
  print('| '.join(var_info + group_info), file = log_stream, end= '\n' if endline else '\r')

def process_one_file(vcf_file, group_map, output_stream, log_stream = sys.stderr):
  upstream = deque(maxlen=args.flank)
  downstream = deque(maxlen=args.flank)
  pysam.set_verbosity(0) # Needed to supress `VariantFile` from pointless warning about an index file
  with pysam.VariantFile(vcf_file) as vcf_handle:
    for record in vcf_handle.fetch():
      
      # Pre-load upstream queue before processing any records
      if len(upstream) < args.flank:
        upstream.append(record)
        continue
      
      # Process one variant and update queues
      current = upstream[0]
      upstream.append(record)
      target_groups, genotypes, primer3_out = process_one_variant(current, group_map, downstream, upstream)
     
      if args.out is not None and args.log is None:
        print_log_stats(log_stream)
      if genotypes is not None: #if a diagnositc variant was found that passed the filter
        for group in target_groups: # It is possible for a variant to be diagnositic for more than one group
          render_output(current, genotypes, reference, primer3_out[group], group, output_stream, downstream, upstream)
      downstream.appendleft(current) # The left side is nearest
      
      
def render_output(record, genotypes, reference, primer3_out, target_group, stream, downstream, upstream):
  # Reformat results
  result = []
  for group in groups: 
    result.append(';'.join([k+str(v) for k, v in genotypes[group].items()]))
    
  # Add reference sequence to the results if reference was supplied
  if reference is not None:
    sequence = extract_seq(record, target_group, downstream, upstream, annotate = args.annotation[0])
    unanot_seq = extract_seq(record, target_group, downstream, upstream, annotate = "minimal")
    if primer3_out is not None and args.annotation[0] != "none":
      pos_to_index_key = flatten([[i]*len(x) for i, x in enumerate(unanot_seq)])
      def annotate(seq, text, pos, left_side = True):
        index = pos_to_index_key[pos]
        item_start_pos = sum([len(x) for x in unanot_seq[:index]])
        item_end_pos = item_start_pos + len(unanot_seq[index]) - 1
        subindex = pos - item_start_pos
        subseq = list(seq[index])
        if left_side:
          if subseq[0] == SNP_DELIM[0]:
            if subindex == 0:
              seq[index] = text + seq[index]
              return
            subindex += 1
          subseq[subindex] = text + subseq[subindex]
        else:
          if subseq[0] == SNP_DELIM[0]:
            if subindex == len(unanot_seq[index]) - 1:
              seq[index] = seq[index] + text
              return
            subindex += 1
          subseq[subindex] = subseq[subindex] + text
        seq[index] = "".join(subseq)
        return
      # annotate right primer
      annotate(sequence, PRIMER_DELIM[0], primer3_out['PRIMER_RIGHT_0'][0] - primer3_out['PRIMER_RIGHT_0'][1] + 1, left_side = True)
      annotate(sequence, PRIMER_DELIM[1], primer3_out['PRIMER_RIGHT_0'][0], left_side = False)
      # annotate crRNA
      spacer_start = get_spacer_start()
      annotate(sequence, CRRNA_DELIM[0], spacer_start, left_side = True)
      annotate(sequence, CRRNA_DELIM[1], spacer_start + args.spacer_len - 1, left_side = False)
      # annotate left primer
      annotate(sequence, PRIMER_DELIM[0], primer3_out['PRIMER_LEFT_0'][0], left_side = True)
      annotate(sequence, PRIMER_DELIM[1], sum(primer3_out['PRIMER_LEFT_0']) - 1, left_side = False)

  if args.format == 'tsv':
    output = [record.chrom, str(record.pos), target_group]
    output.extend(result)
    if reference is not None:
      spacer_start = get_spacer_start()
      spacer_seq = "".join(unanot_seq)[spacer_start:spacer_start + args.spacer_len]
      output.append("".join(sequence))
      output.append(spacer_seq)
    if primer3_out is not None:
      # primer3_data = [str(round(primer3_out[n], 3)) if type(primer3_out[n]) == float  else str(primer3_out[n]) for n in primer3_col_names]
      primer3_data = [str(primer3_out[n]) for n in primer3_col_names]
      output.extend(primer3_data)
    print('\t'.join(output), file = stream)
  elif args.format == 'csv':
    raise Exception('Not yet implemented')
  elif args.format == 'vcf':
    raise Exception('Not yet implemented')
  elif args.format == 'human':
    raise Exception('Not yet implemented')
  else:
    raise Exception('Not yet implemented')


def render_header(stream):
  if args.format in ('tsv', 'csv'):
    if args.format == 'tsv':
      out_sep = '\t'
    else:
      out_sep = ','
    if args.reference is None:
      print("chormosome", "position", "target", *groups, sep=out_sep, file = stream)
    elif args.primer3:
      print("chormosome", "position", "target", *groups, "reference", "spacer", *[primer3_col_key[n] for n in primer3_col_names], sep=out_sep, file = stream)
    else:
      print("chormosome", "position", "target", *groups, "reference", "spacer", sep=out_sep, file = stream)
  elif args.format == 'vcf':
    raise Exception('Not yet implemented')
  elif args.format == 'human':
    raise Exception('Not yet implemented')
  else:
    raise Exception('Not yet implemented')


#https://stackoverflow.com/questions/22264504/using-python-with-statement-with-sys-stdout
@contextmanager
def writer(file_path = None, default_stream = sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "w")
    yield file_handle
    if file_path != None:
      file_handle.close()


# Entry point into the program and major steps
if __name__ == "__main__":
  
  # Parse arguments
  args = parser.parse_args()
  group_map = get_group_map(args.metadata)
  group_counts = Counter(group_map.values())
  reference = get_reference(args.reference)
  if args.groups is None:
    groups = list({k for k, v in group_counts.items() if v >= args.min_samples})
  else:
    groups = args.groups
  
  # Validate arguments
  if args.primer3 and reference is None:
    parser.error('A reference is required (--reference option) to pick primers.')
  for group, count in group_counts.items():
    if group in groups:
      if count < args.min_samples:
        parser.error(f"The group '{group}' has fewer samples ({count}) than --min_samples ({args.min_samples}), so no sites would be returned.")
  for group in groups:
    if group not in group_map.values():
      parser.error(f"The group '{group}' is not present in the metadata.")
    
  # Set up stat counters
  variant_stat_names = ['Processed', 'Undiagnostic', 'Missing data', 'Low map qual', 'Unconserved', 'No primers']
  variant_stats = {k: 0 for k in variant_stat_names}
  group_counts = {k: 0 for k in groups}
  
  # Set up primer3 columns to print
  primer3_col_names = [
    'PRIMER_PAIR_0_PRODUCT_SIZE',
    'PRIMER_PAIR_0_PENALTY',
    'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE',
    'PRIMER_LEFT_0_PENALTY', 'PRIMER_RIGHT_0_PENALTY',
    'PRIMER_LEFT_0_TM','PRIMER_RIGHT_0_TM',
    'PRIMER_LEFT_0_GC_PERCENT','PRIMER_RIGHT_0_GC_PERCENT',
    'PRIMER_LEFT_0_SELF_ANY_TH', 'PRIMER_RIGHT_0_SELF_ANY_TH',
    'PRIMER_LEFT_0_SELF_END_TH', 'PRIMER_RIGHT_0_SELF_END_TH',
    'PRIMER_LEFT_0_HAIRPIN_TH', 'PRIMER_RIGHT_0_HAIRPIN_TH',
    'PRIMER_LEFT_0_END_STABILITY', 'PRIMER_RIGHT_0_END_STABILITY',
    'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH'
  ]
  primer3_col_key = {n: n.replace("PRIMER_", "").replace("_0", "").lower() for n in primer3_col_names}
  
  with writer(args.out, sys.stdout) as output_stream, writer(args.log, sys.stderr) as log_stream:
    # print header if needed
    render_header(output_stream)
    if args.out is not None and args.log is None:
      print_log_stats_header(log_stream)
    # Process each input file
    for vcf_file in args.vcfs:
      process_one_file(vcf_file, group_map, output_stream, log_stream)
    if args.out is None or args.log is not None:
      print_log_stats_header(log_stream)
    print_log_stats(log_stream, endline = True)

