import gzip
import multiprocessing
import shutil
import os
import sys
from .shared import *
from contextlib import contextmanager
import pdb

from pudb.remote import set_trace


@contextmanager
def stream_writer(file_path=None, default_stream=sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "a")
    yield file_handle
    if file_path is not None and file_handle is not None:
        file_handle.close()


def safe_print(alignments, tsv_rows, out_align, out_tsv, counter, lock):
    """ Function to safely print alignments to the terminal in parallel """
    # Acquire lock
    with lock:
        # Iterate through alignments and print, increment counter
        if out_align is not None and out_tsv is not None:
            for align, row in zip(alignments, tsv_rows):
                print(align, file=out_align, flush=True)
                print(row, file=out_tsv, flush=True)
                counter.value += 1
        elif out_tsv is not None:
            for row in tsv_rows:
                print(row, file=out_tsv, flush=True)
                counter.value += 1
        elif out_align is not None:
            for align in alignments:
                print(align, file=out_align, flush=True)
                counter.value += 1


def render_output_part(kmerfile, out_align, out_tsv, counter, lock,
                       start=None, end=None, print_block=100, ingroup=None, find_primers=False):
    """ Helper function to render the contents of a part of the kmer file """
    # Get alignment stream
    alignments = alignmentStream(kmerfile, start, end, ingroup)
    alignments_to_print = []
    row_to_print = []
    print_block_counter = 0
    with stream_writer(out_align, None) as out_align_stream,\
         stream_writer(out_tsv, sys.stdout) as out_tsv_stream:
        for alignment in alignments:
            # Run primer3 (NOTE: ideally this would run on its own step in krisp.main, not here)
            if find_primers:
                if not alignment.find_primers():
                    continue
            # Add to alignments list
            if out_align_stream is not None:
                alignments_to_print.append(alignment.render_alignment())
            # Add TSV rows to list
            if out_tsv_stream is not None:
                row_to_print.append(alignment.render_tsv())
            #set_trace(term_size=(80, 60))
            # Print
            print_block_counter += 1
            if print_block_counter >= print_block:
                # Print alignments
                safe_print(alignments_to_print, row_to_print, out_align_stream, out_tsv_stream, counter, lock)
                alignments_to_print = []
                row_to_print = []
                print_block_counter = 0
        if print_block_counter > 0:
            safe_print(alignments_to_print, row_to_print, out_align_stream, out_tsv_stream, counter, lock)


def render_output(kmerfile, out_align=None, out_tsv=sys.stdout, cores=1, print_block=10, ingroup=None, find_primers=False):
    """ Function to print alignments from a file in parallel

    Parameters
    ----------

    kmerfile : str
        The input kmerfile to render

    out_align : str or stream
        The path to the alignment output file or stream to write to (e.g. sys.stdout)

    out_tsv : str or stream
        The path to the TSV output file or stream to write to (e.g. sys.stdout)

    cores : int
        The number of processors to use

    print_block : int
        The number of results to buffer in memory before attempting to acquire the lock and write the results.

    ingroup : list of str
        The names of the in-group samples

    Returns
    -------
    int
        The number of alignments written

    """

    # Get the start and end iteration points for kmerfile
    fptr_start_end = list(splitFilePtrs(kmerfile, num_sections=cores))

    # Create a counter and lock
    counter = multiprocessing.Value('i', 0)
    lock = multiprocessing.Lock()

    # Create a list of job args
    job_args = []
    for start, end in fptr_start_end:
        job_args.append((kmerfile, out_align, out_tsv, counter, lock, start, end, print_block, ingroup, find_primers))

    # Start a group of processes
    processes = []
    for job in job_args:
        p = multiprocessing.Process(target=render_output_part,
                                    args=job)
        p.start()
        processes.append(p)

    # Call join on processes
    for p in processes:
        p.join()
    return counter.value



# def writeAlignments(kmerfile, output, pid, start=None, end=None, ingroup=None):
#     """ Helper function to write alignments from kmer file """
#     # Get alignment stream
#     alignments = alignmentStream(kmerfile, start, end, ingroup)
#
#     # Open for writing and write output
#     with gzip.open(output, 'wt') as fout:
#         found = 0
#         for i, alignment in enumerate(alignments):
#             # Print allignment
#             fout.write(f"> process_{pid}_{i}\n")
#             fout.write(f"{str(alignment)}\n")
#             found += 1
#         fout.flush()
#     return found
#
#
# def writeAlignmentsParallel(kmerfile, output, parallel, workdir=None, ingroup=None):
#     """ Function to write alignments from a kmer file in parallel
#
#     Parameters
#     ----------
#
#     kmerfile : str
#         The input kmerfile to convert to an alignment file
#
#     output : str
#         The output alignment file
#
#     parallel : int
#         The number of processors to use
#
#     workdir : str (default None)
#         The working directory to place temporary files
#
#     Returns
#     -------
#     int
#         The number of alignments written
#
#     """
#
#     # Get the start and end iteration points for kmerfile
#     fptr_start_end = list(splitFilePtrs(kmerfile, parallel))
#
#     # Create a list of temporary output files
#     tmp_files = [tmpFile(workdir) for i in range(parallel)]
#
#     # Create a list of job args
#     job_args = []
#     for i in range(parallel):
#         start, end = fptr_start_end[i]
#         pid = i
#         tmp = tmp_files[i]
#         job_args.append((kmerfile, tmp, pid, start, end, ingroup))
#
#     # Run processes in pool
#     with multiprocessing.Pool(parallel) as pool:
#         counts = pool.starmap(writeAlignments, job_args)
#
#     # Merge tmpfiles into output
#     if output.endswith('.gz'):
#         with open(output, "wb") as fout:
#             for tmp in tmp_files:
#                 with open(tmp, "rb") as fin:
#                     shutil.copyfileobj(fin, fout)
#     else:
#         with open(output, "w") as fout:
#             for tmp in tmp_files:
#                 with gzip.open(tmp, "rt") as fin:
#                     shutil.copyfileobj(fin, fout)
#     return sum(counts)
