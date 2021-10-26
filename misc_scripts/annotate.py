#! /bin/env python3
from krisp.shared import *
import sys
import argparse
from grepFastaSequence import findInFasta


def seqToDotAlignment(ref, seq):
    result = []
    for r, s in zip(ref, seq):
        if r == s:
            result.append('.')
        else:
            result.append(s)
    return ''.join(result)


def dotToSeqAlignment(ref, seq):
    result = []
    for r, s in zip(ref, seq):
        if s == '.':
            result.append(r)
        else:
            result.append(s)
    return ''.join(result)    


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Annotate alignments from krisp output",
            epilog="Hope this helped. GET BACK TO WORK!",
            prog="annotate")
    parser.add_argument(
            "alignment",
            type=str,
            help="Alignment file to annotate")
    parser.add_argument(
            "-r",
            "--ref",
            nargs="+",
            type=str,
            help="Reference fasta files used to generate the alignment")
    parser.add_argument(
            "-p",
            "--padding",
            type=int,
            default=0,
            help="Add upstream and downstream padded region for each sequence")
    parser.add_argument(
            "-o",
            "--output",
            type=str,
            help="Write results to a file")
    args = parser.parse_args(sys.argv[1:])

    # Set output
    if args.output is None:
        args.output = sys.stdout
    else:
        args.output = open(args.output, "w")

    # Get mapping of basename to ref_file
    base_to_file = {simplename(f) : f for f in args.ref}

    # Iterate through each alignment
    ref_seq = None
    ref_buffered_seq = None
    count = 0
    for line in open(args.alignment):
        # Check if beginning of new alignment
        line = line.strip()
        if line.startswith('>'):
            ref_seq = None
            ref_buffered_seq = None
            count += 1
            print(f"# Alignment {count}", file=args.output)
            continue
        elif line.startswith('{'):
            continue

        # Found first sequence of an alignment
        seq, tags = line.split(':')
        seq = seq.strip()

        # Remove dots
        if ref_seq is None:
            ref_seq = seq
        else:
            seq = dotToSeqAlignment(ref_seq, seq)

        # Parse tags
        tags = tags.split(';')
        tags = [tag.split('(')[0].strip() for tag in tags]

        # get seqs, lables
        first = True
        for tag in tags:
            filename = base_to_file[tag]
            for name, start, end, sequence in findInFasta(filename,
                                                    seq,
                                                    args.padding):
                # Make tag
                short_name = filename.split('/')[-1]
                is_reverse = sequence.islower()
                strand = "forward" if not is_reverse else "reverse" 
                label = f"{short_name} | chrom={name} | start={start} | end={end} | strand={strand}"
                sequence = sequence.upper()

                #print(f"> {label}")
                if ref_buffered_seq is None:
                    print(f"{sequence} | {label}", file=args.output)
                    ref_buffered_seq = sequence
                else:
                    dot = seqToDotAlignment(ref_buffered_seq, sequence)
                    print(f"{dot} | {label}", file=args.output)

    # Close file
    if args.output is not sys.stdout:
        args.output.close()