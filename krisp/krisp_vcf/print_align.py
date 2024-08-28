import shutil
import math
from collections import defaultdict
from prettytable import PrettyTable

class Annotation:

    def __init__(self, name, seq, start):
        self.name = name
        self.seq = seq
        self.start = start


def cumulative(lists):
    """
    https://www.geeksforgeeks.org/python-program-to-find-cumulative-sum-of-a-list/
    """
    length = len(lists)
    cu_list = [sum(lists[0:x:1]) for x in range(0, length + 1)]
    return cu_list[1:]


def _mask_same(seqs, ref):
    """If loci are the same, letter will be replaced with a period

    Parameters
    ----------
    seqs : dict of list of str
        sequences to align, named by group
    ref : list of str
        the reference used to call the variants

    Returns
    -------
    list of str
        Modified version of seqs with equivalent elements replaced with "."

    """
    for group in seqs.keys():
        for row_index in range(len(seqs[group])):
            if seqs[group][row_index] == ref[row_index]:
                middle = len(ref[row_index])
                seqs[group][row_index] = "." * middle

    return seqs


def pos_to_chunk_index(pos, ref):
    pos_to_chunk_key = {p - 1: i for i, p in enumerate(cumulative([len(c.strip()) for c in ref]))}
    for p, i in pos_to_chunk_key.items():
        if p >= pos:
            return {'chunk': i, 'offset': abs(p - pos)}


def _pad_sequences(seqs, ref, annots):
    """Adds padding for all sequences, with spacing or dashes
    Parameters
    ----------
    seqs : dict of list of str
        sequences to align, named by group
    ref : list of str
        list of alleles in the reference genome

    Return
    ------
    list of str
        Modified ref where indels are labelled "-", up to max. column length

    dict of list of str
        Modified seqs where indels are labelled "-", up to max. column length,
        any ambiguous or heterozygous loci are padded with " " to max. length

     """

    def _pad_ref_and_seqs(col_i, width, pad_str):
        ref[col_i] = ref[col_i].center(width, pad_str)
        for seq_name in seqs.keys():
            if seqs[seq_name][col_i] == " ": #used for printing oligos sequences
                seqs[seq_name][col_i] = seqs[seq_name][col_i].center(width, " ")
            else:
                seqs[seq_name][col_i] = seqs[seq_name][col_i].center(width,
                                                                     pad_str)
        return None

    # Pad reference and sequences to equal column widths
    for col_index in range(len(ref)):
        column = [ref[col_index]] + [seq[col_index] for seq in seqs.values()]
        max_width = max([len(x) for x in column])
        if any([">" in x or "/" in x for x in column]):
            _pad_ref_and_seqs(col_index, max_width, pad_str=" ")
        else:
            _pad_ref_and_seqs(col_index, max_width, pad_str="-")

    # Split annotation seq names into columns for printing later
    col_widths = [len(c) for c in ref]
    annot_out = [" " * len(c) for c in ref]
    for annot in annots:
        start = pos_to_chunk_index(annot.start, ref)
        end = pos_to_chunk_index(annot.start + len(annot.seq) - 1, ref)
        annot_col_widths = [col_widths[i] for i in range(start['chunk'], end['chunk'] + 1)]
        print_len = sum(annot_col_widths)
        name = ' ' + annot.name + ' '
        text = '└' + name.center(print_len - 2, '─') + '┘'
        text_inter = iter(text)
        text_cols = ["".join(next(text_inter) for i in range(w)) for w in annot_col_widths]
        for text_index, ref_index in enumerate(range(start['chunk'], end['chunk'] + 1)):
            annot_out[ref_index] = text_cols[text_index]

    return seqs, ref, annot_out


def _print_align(seqs, ref, annot_text, groups, ref_name="Reference"):
    """Adjusts labels, prints reference and sequences
    """


    def _print_one_line(seqs, ref, groups, ref_name="Reference"):
        """Repeated function to print the output format"""
        # Get counts of each group to print
        group_counts = {g: str(len(v)) for g, v in groups.items()}
        max_count_len = max([len(n) for n in group_counts.values()])

        # Get max length of labels
        labels = {k: f'{k} ({group_counts[k]})' if k in group_counts else f'{k}' for k in seqs.keys()}
        max_len = max([len(label) for label in list(labels.values()) + [ref_name]])

        # Print alignment
        output = []
        ref_name = ref_name.rjust(max_len, " ")
        #import pdb; pdb.set_trace()


        output.append(f"{ref_name}: " + ''.join(ref))
        for group_name, seq in seqs.items():
            group_text = labels[group_name].rjust(max_len, " ")
            output.append(f"{group_text}: " + ''.join(seq))

        # Print annotation names
        output.append(' ' * (max_len + 2) + ''.join(annot_text))
        # ref_len = sum([len(x) for x in ref])
        # annot_line = ([" "] * ref_len)
        # for annot in annots:
        #     size =
        #     for index, nucleotide in enumerate(annot.seq):
        #         annot_line[annot.start + index] = nucleotide

        return output

    def split(x, f):
        """
        https://stackoverflow.com/questions/19597097/python-equivalent-of-r-split-function
        """
        res = defaultdict(list)
        for v, k in zip(x, f):
            res[k].append(v)
        return res
#   align_width: consistent label length so sequences line up
#   col_widths: find length of columns in ref
#   row_index: length of character divided by aligned width, rounded down
#   chunked_ref: sort row_index into the assigned rounded groups
#   chunked_seqs: same as chunked_ref, but for dict of list (dict comprehension)
#   for loop makes the columns aesthetic
    term_width = shutil.get_terminal_size().columns
    labels = list(seqs.keys()) + [ref_name]
    label_width = max([len(seq) for seq in labels])
    align_width = term_width - label_width - 5
    col_widths = [len(x) for x in ref]
    row_index = [math.floor(x/align_width) for x in cumulative(col_widths)]
    chunked_ref = split(ref, row_index)
    chunked_seqs = {k: split(v, row_index) for k, v in seqs.items()}
    output = []
    for index in range(len(chunked_ref)):
        row_seqs = {k: v[index] for k, v in chunked_seqs.items()}
        output.extend(_print_one_line(row_seqs, chunked_ref[index], groups, ref_name=ref_name))
    return output


def format_seq_annot(annots, ref):
    """Go through annotated sequences and overwrite to add primer data

    Parameters:
    -----------
    annots: list of Annotation
        primers identified by Primer3/ package
    ref: list of str
        list of alleles in the reference genome
    Return:
         list(?)
    """
    ref_len = sum([len(x) for x in ref])
    output = ([" "] * ref_len)
    for annot in annots:
        start = pos_to_chunk_index(annot.start, ref)
        for index, nucleotide in enumerate(annot.seq):
            output[start['chunk'] + index] = nucleotide
    return output


def _render_primer3_stats(p3):
    left_stats = {k[14:]: v for k, v in p3.items() if 'PRIMER_LEFT_0_' in k}
    right_stats = {k[15:]: v for k, v in p3.items() if 'PRIMER_RIGHT_0_' in k}
    pair_stats = {k[14:]: v for k, v in p3.items() if 'PRIMER_PAIR_0_' in k}

    def render_col_names(names):
        return [x.title().replace('_', ' ') for x in names]

    def render_col_values(names):
        return [str(round(x, 5)) if isinstance(x, float) else x for x in names]

    primer_table = PrettyTable(['Direction'] + render_col_names(left_stats.keys()))
    primer_table.add_row(['Forward'] + render_col_values(left_stats.values()))
    primer_table.add_row(['Reverse'] + render_col_values(right_stats.values()))
    primer_table.align = 'l'

    pair_table = PrettyTable(render_col_names(pair_stats.keys()))
    pair_table.add_row(render_col_values(pair_stats.values()))
    pair_table.align = 'l'

    output = '\nPrimer statistics:\n' + \
             primer_table.get_string(border=False) + \
             '\n\nPair statistics:\n' + \
             pair_table.get_string(border=False)
    return output


def render_variant(seqs, ref, p3, groups, annots=None):
    """Displaying diagnostic variant in human-readable form

    Parameters
    ----------
    seqs : dict of list of str
        sequences to align, named by group
    ref : list of str
        the reference used to call the variants
    p3 : dict
        raw primer3 output
    annots : list of Annotation
        sequences to annotate alignment, keys are sequences, values are start
        position in the alignment

    Returns
    -------
    list of str
        one str for each line to print
    """

    # Make alignment
    seqs = _mask_same(seqs, ref)
    if annots is not None:
        seqs["oligos"] = format_seq_annot(annots, ref)
    seqs, ref, annot_text = _pad_sequences(seqs, ref, annots)
    output = _print_align(seqs, ref, annot_text, groups)

    # Make Primer3 output
    output += [_render_primer3_stats(p3)]

    return output


if __name__ == "__main__":
    ref = ["A", "C", "GTC", "T", "T", "G", "G", "C", "C", "A", "C", "GTC", "T",
           "T", "G", "G", "C", "C", "A", "C", "GTC", "T", "T", "G", "G", "C",
           "C", "A", "C", "GTC", "T", "T", "G", "G", "C", "C", "A", "C", "GTC",
           "T", "T", "G", "G", "C", "C", "A", "C", "GTC", "T", "T", "G", "G",
           "C", "C", "A", "C", "GTC", "T", "T", "G", "G", "C", "C"]
    seqs = {"seq 1": ["A", "G", "GTC", "T/G", "T", "GC", "G", "A", "C", "A",
                      "G", "GTC", "T/G", "T", "GC", "G", "A", "C", "A", "G",
                      "GTC", "T/G", "T", "GC", "G", "A", "C", "A", "G", "GTC",
                      "T/G", "T", "GC", "G", "A", "C", "A", "G", "GTC", "T/G",
                      "T", "GC", "G", "A", "C", "A", "G", "GTC", "T/G", "T",
                      "GC", "G", "A", "C", "A", "G", "GTC", "T/G", "T", "GC",
                      "G", "A", "C"],
            "seq 2": ["A", "G", "G", "T", "", "G", "<123G12?>", "A", "C", "A",
                      "G", "G", "T", "", "G", "<123G12?>", "A", "C", "A", "G",
                      "G", "T", "", "G", "<123G12?>", "A", "C", "A", "G", "G",
                      "T", "", "G", "<123G12?>", "A", "C", "A", "G", "G", "T",
                      "", "G", "<123G12?>", "A", "C", "A", "G", "G", "T", "",
                      "G", "<123G12?>", "A", "C", "A", "G", "G", "T", "", "G",
                      "<123G12?>", "A", "C"]}
    annotation = [Annotation("Left primer", "ATTGCATGA", 33),
                  Annotation("crRNA", "TGGTCCATGAT", 45),
                  Annotation("Right primer", "ATTGTAAC", 12)]
    render_variant(seqs, ref, annots=annotation)


