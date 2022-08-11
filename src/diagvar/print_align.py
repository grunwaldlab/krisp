import shutil
import math
from collections import defaultdict


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


def _pad_sequences(seqs, ref):
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

    for col_index in range(len(ref)):
        column = [ref[col_index]] + [seq[col_index] for seq in seqs.values()]
        max_width = max([len(x) for x in column])
        if any([">" in x or "/" in x for x in column]):
            _pad_ref_and_seqs(col_index, max_width, pad_str=" ")
        else:
            _pad_ref_and_seqs(col_index, max_width, pad_str="-")
    return seqs, ref


def _print_align(seqs, ref, ref_name="Reference"):
    """Adjusts labels, prints reference and sequences
    """

    def cumulative(lists):
        """
        https://www.geeksforgeeks.org/python-program-to-find-cumulative-sum-of-a-list/
        """
        length = len(lists)
        cu_list = [sum(lists[0:x:1]) for x in range(0, length + 1)]
        return cu_list[1:]

    def _print_one_line(seqs, ref, ref_name="Reference"):
        """Repeated function to print the output format"""
        labels = list(seqs.keys()) + [ref_name]
        max_len = max([len(seq) for seq in labels])
        ref_name = ref_name.rjust(max_len, " ")
        print(f"{ref_name} : " + ''.join(ref))
        for group_name, seq in seqs.items():
            group_name = group_name.rjust(max_len, " ")
            print(f"{group_name} : " + ''.join(seq))

    def split(x, f):
        """
        https://stackoverflow.com/questions/19597097/python-equivalent-of-r-split-function
        """
        res = defaultdict(list)
        for v, k in zip(x, f):
            res[k].append(v)
        return res
#   term_width: get the width of the terminal
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
    for index in range(len(chunked_ref)):
        row_seqs = {k: v[index] for k, v in chunked_seqs.items()}
        _print_one_line(row_seqs, chunked_ref[index], ref_name=ref_name)
        print()


def format_seq_annot(primer, ref):
    """Go through annotated sequences and overwrite to add primer data

    Parameters:
    -----------
    primer: dict of list of str
        primers identified by Primer3/ package
    ref: list of str
        list of alleles in the reference genome
    Return:
         list(?)
    """
    output = ([" "] * len(ref))
    for seq, start in primer.items():
        for index, nucleotide in enumerate(seq):
            output[start + index] = nucleotide

    return output


def render_variant(seqs, ref, annots=None):
    """Displaying diagnostic variant in human-readable form

    Parameters
    ----------
    seqs : dict of list of str
        sequences to align, named by group
    ref : list of str
        the reference used to call the variants
    annots : dict of str, int
        sequences to annotate alignment, keys are sequences, values are start
        position in the alignment
    """
    seqs = _mask_same(seqs, ref)
    if annots is not None:
        seqs["oligos"] = format_seq_annot(annots, ref)
    seqs, ref = _pad_sequences(seqs, ref)
    _print_align(seqs, ref)

    return None


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
    annotation = {'ATTGCATGA': 33, 'TGGTCCATGAT': 45, 'ATTGTAAC': 12}
    #TODO: add information about forward and reverse primers
    render_variant(seqs, ref, annots=annotation)
