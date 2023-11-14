#! /bin/env python3
import sys
import argparse
import time
import multiprocessing
import subprocess
from ..kstream import kstream
from .intersectAmplicons import mergeFiles
from .outputAlignments import render_output
from .shared import *
from colorama import Fore, Back, Style
from .filterAlignments import filterAlignments
from .Amplicon import ConservedEndAmplicons


def extractSortedKmers(fasta, primer_left, primer_right, ampl_len, output,
                       sortmem, parallel=1, verbose=True, omit=True):
    """  Fastafile -> kmers written to output file """
    # Get kmers from the fasta file
    if omit:
        kmers = kstream(fasta,
                        kmers=ampl_len,
                        disallow="Nn",
                        complements=True,
                        omitsoft=True,
                        split=[primer_left, -primer_right],
                        sort=True,
                        sortmem=sortmem,
                        sortcols=[0, 2],
                        sortnp=parallel,
                        parallel=parallel)
    else:
        kmers = kstream(fasta,
                        kmers=ampl_len,
                        disallow="Nn",
                        complements=True,
                        mapsoft=True,
                        split=[primer_left, -primer_right],
                        sort=True,
                        sortmem=sortmem,
                        sortcols=[0, 2],
                        sortnp=parallel,
                        parallel=parallel)


    # Write kmers to output
    if verbose:
        # Get start time
        start_t = time.time()

        # Print start message
        message = (f"Extracting {ampl_len}-mers from {fasta} "
                   f"and saving to {output}")
        print(message, end='\n', file=sys.stderr)

        # Write kmers and get count
        found = kmers.write(output)

        # Print end message
        end_t = time.time()
        end_message = (f"=> Extracted and sorted {found:,} {ampl_len}-kmers"
                       f" from {fasta} in {prettyTime(end_t-start_t)}")
        print(Fore.GREEN + end_message + Style.RESET_ALL, file=sys.stderr)
    else:
        # Write kmers
        kmers.write(output)


def sortedKmerJob(input_queue):
    """ Process job which calls kmerjob on args in queue """
    args = input_queue.get()
    while args is not None:
        # run job
        extractSortedKmers(*args)
        args = input_queue.get()


def sortedKmersSerial(files, outputs, ampl_len, primer_left, primer_right,
                      verbose=True, omit=True):
   for filename, outfile in zip(files, outputs):
        # Call base function to extract and sort args
        extractSortedKmers(filename, primer_left, primer_right, ampl_len,
                           outfile, "80%", 1, verbose, omit)


def sortedKmersParallel(files, outputs, ampl_len, primer_left, primer_right,
                        parallel=1, verbose=True, omit=True):
    """ Coverts a batch of files into sorted kmer files """
    # Determine the number of cores to give to each job
    job_cores = []
    jobs = len(files)
    if parallel < jobs:
        job_cores = [1] * jobs
    else:
        for i, f in enumerate(files):
            np = (parallel // jobs) + (i < parallel % jobs)
            job_cores.append(np)

    # First add jobs to a multiprocessing queue
    job_queue = multiprocessing.Queue()
    sortmem = f"{80//len(job_cores)}%"
    for filename, outfile, cores in zip(files, outputs, job_cores):
        # Get args for finding kmers
        args = (filename, primer_left, primer_right, ampl_len,
                outfile, sortmem, cores, verbose, omit)
        # Add args to job queue
        job_queue.put(args)
    for i in range(jobs):
        # Add end signal to job queue
        job_queue.put(None)

    # Spawn jobs
    processes = []
    for i in range(min([parallel, jobs])):
        # Spawn process
        p = multiprocessing.Process(target=sortedKmerJob, args=(job_queue,))
        # Start process and add to process queue
        p.start()
        processes.append(p)

    # Wait for jobs to finish
    for p in processes:
        p.join()


def main():
    """ Parse command line args """
    parser = argparse.ArgumentParser(
            description="Find diagnostic alignments for a set of fasta files",
            prog="krisp",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("files", nargs="+", type=str, metavar='PATH',
                        help="Fasta file to read. .gz, .bz2")
    parser.add_argument("--outgroup", nargs="*", type=str, default=[], metavar='PATH',
                        help="Outgroup Fasta files. To be amplified, but not detected")
    parser.add_argument("-c", "--conserved", type=int, metavar='INT',
                        help="Length of conserved regions on ends of amplicon")
    parser.add_argument("--conserved-left", type=int, metavar='INT',
                        help="Length of conserved region on left of amplicon")
    parser.add_argument("--conserved-right", type=int, metavar='INT',
                        help="Length of conserved region on right of amplicon")
    parser.add_argument("-d", "--diagnostic", type=int, metavar='INT',
                        help="Diagnostic region length for amplicon")
    parser.add_argument("-a", "--amplicon", type=int, metavar='INT',
                        help="Total amplicon length")
    parser.add_argument("--omit-soft", action="store_true",
                        help="Omit softmasked nucleotides")
    parser.add_argument("--cores", type=int, default=1, metavar='INT',
                        help="Total number of processors to utilize. (default: %(default)s)")
    parser.add_argument("--dot-alignment", action="store_true",
                        help="Output as dot-based alignments")
    parser.add_argument("-o", "--out_align", type=str, metavar='PATH',
                        help="Write results as human-readable alignments to a file (gzip supported). (default: do not write alignment output)")
    parser.add_argument("-s", "--out_csv", type=str, metavar='PATH',
                        help="Write results to as a CSV (comma-separated value) file (gzip supported). (default: print to screen (stdout))")
    parser.add_argument("-w", "--workdir", type=str, metavar='PATH',
                        help="Work directory to place temporary files")
    parser.add_argument("-p", "--primer3", action=argparse.BooleanOptionalAction,
                        help="Work directory to place temporary files")
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
    parser.add_argument("--verbose", action="store_true",
                        help="Print runtime information to sys.stderr")
    args = parser.parse_args(sys.argv[1:])

    # Set conserved left, conserved right and amplicon
    if args.amplicon is not None:
        if args.diagnostic is not None:
            # Determine conserved region lengths
            args.conserved = (args.amplicon - args.diagnostic) // 2
            args.conserved_left = args.conserved
            args.conserved_right = args.conserved
        elif args.conserved is not None:
            # Determine diagnostic region length
            args.diagnostic = args.amplicon - 2 * args.conserved
            args.conserved_left = args.conserved
            args.conserved_right = args.conserved
        elif (args.conserved_left is not None) and (args.conserved_right is not None):
            # Set diagnostic length
            args.diagnostic = (args.amplicon - args.conserved_left - args.conserved_right)
        else:
            print("ERROR: Could not deduce input parameters", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)
    elif args.diagnostic is not None:
        if args.conserved is not None:
            # Determine amplicon length
            args.amplicon = args.diagnostic + 2 * args.conserved
            args.conserved_left = args.conserved
            args.conserved_right = args.conserved
        elif (args.conserved_left is not None) and (args.conserved_right is not None):
            # Set diagnostic length
            args.amplicon = args.diagnostic + args.conserved_left + args.conserved_right
        else:
            print("ERROR: Could not deduce input parameters", file=sys.stderr)
            parser.print_help(sys.stderr)
            sys.exit(1)
    else:
        print("ERROR: Could not deduce input parameters", file=sys.stderr)
        parser.print_help(sys.stderr)
        sys.exit(1)            

    # Set output format
    ConservedEndAmplicons.ENABLE_DOT = args.dot_alignment

    # Set Primer3 parameters (NOTE: should be passed a different way, but this was easy)
    primer3_arg_names = ('tm', 'gc', 'primer_size', 'amp_size', 'max_sec_tm',
                         'min_bases', 'gc_clamp', 'max_end_gc')
    ConservedEndAmplicons.P3_ARGS = {k: v for k, v in vars(args).items() if k in primer3_arg_names}

    # Create a temporary directory to work in
    with tempfile.TemporaryDirectory(dir=args.workdir) as tmpdir:
        # Time if verbose
        if args.verbose:
            start_t = time.time()
            print("Finding kmer-based diagnostic regions for:",
                  file=sys.stderr)
            for i, filename in enumerate(args.files):
                print(f"({i}) {filename}", file=sys.stderr)
            print("With this as an outgroup:", file=sys.stderr)
            for i, filename in enumerate(args.outgroup):
                print(f"({i}) {filename}", file=sys.stderr)
            print(file=sys.stderr)

        # For every file, extract kmers and save to an individual file
        input_files = args.files + args.outgroup
        kmer_files = []
        for filename in input_files:
            kmer_name = (f"{tmpdir}/{basename(Path(filename).name)}"
                         f".{args.amplicon}mers")
            kmer_files.append(kmer_name)

        # Get sorted kmers
        if args.cores > 1:
            sortedKmersParallel(input_files, kmer_files, args.amplicon,
                                args.conserved_left, args.conserved_right, args.cores,
                                verbose=args.verbose, omit=args.omit_soft)
        else:
            sortedKmersSerial(input_files, kmer_files, args.amplicon,
                                args.conserved_left, args.conserved_right,
                                verbose=args.verbose, omit=args.omit_soft)

        # Merge kmer files into a single file
        result = f"{tmpdir}/merged_file.txt"
        mergeFiles(kmer_files, result, args.cores, tmpdir, args.verbose)

        # Print start of alignment building
        if args.verbose:
            print(file=sys.stderr)
            print("Filtering for diagnostic regions ... ", end='\n', file=sys.stderr)

        # Find diagnostic sets matching this pattern
        if (args.amplicon > args.conserved_left + args.conserved_right):
            # Create a frozen set from ingroup names
            ingroup_set = frozenset([simplename(f) for f in args.files])
            # Get filtered alignments
            filtered_result = f"{tmpdir}/filtered.txt"
            filterAlignments(result, filtered_result, ingroup_set)
            # Set filtered result as result
            result = filtered_result

        # Print start of alignment building
        if args.verbose:
            print(file=sys.stderr)
            print("Rendering output ... ", end='\n', file=sys.stderr)

        # Print alignments
        found = None
        ingroup = None
        if len(args.outgroup):
            ingroup = [simplename(f) for f in args.files]
        found = render_output(result,
                              out_align=args.out_align,
                              out_csv=args.out_csv,
                              cores=args.cores,
                              print_block=1000,
                              ingroup=ingroup,
                              find_primers=args.primer3)

        # Print end message
        if args.verbose:
            # Get end time
            end_t = time.time()
            end_message = (f"=> Found {found:,} regions"
                           f" in {prettyTime(end_t-start_t)}")
            print(Fore.GREEN + end_message + Style.RESET_ALL, file=sys.stderr)


if __name__ == "__main__":
    main()
