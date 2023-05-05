import gzip
import multiprocessing
import shutil
import os
import sys
from .shared import *
from contextlib import contextmanager


@contextmanager
def stream_writer(file_path=None, default_stream=sys.stdout):
    file_handle = default_stream if file_path is None else open(file_path, "w")
    yield file_handle
    if file_path is not None:
        file_handle.close()


def safe_print(alignments, output_stream, counter, lock):
    """ Function to safely print alignments to the terminal in parallel """
    # Acquire lock
    with lock:
        # Iterate through alignments and print, increment counter
        for align in alignments:
            print(align, file=output_stream, flush=True)
            counter.value += 1


def render_output_part(kmerfile, output, counter, lock,
                       start=None, end=None, print_block=10, ingroup=None):
    """ Helper function to render the contents of a part of the kmer file """
    # Get alignment stream
    alignments = alignmentStream(kmerfile, start, end, ingroup)
    alignments_to_print = []
    with stream_writer(output, sys.stdout) as output_stream:
        for alignment in alignments:
            # Add to alignments list
            alignments_to_print.append(str(alignment))
            # Print
            if len(alignments_to_print) >= print_block:
                # Print alignments
                safe_print(alignments_to_print, output_stream, counter, lock)
                alignments_to_print = []
        if len(alignments_to_print):
            safe_print(alignments_to_print, output_stream, counter, lock)


def render_output(kmerfile, output=sys.stdout, cores=1, print_block=10, ingroup=None):
    """ Function to print alignments from a file in parallel

    Parameters
    ----------

    kmerfile : str
        The input kmerfile to render

    output : str or stream
        The path to the output file or stream to write to (e.g. sys.stdout)

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
        job_args.append((kmerfile, output, counter, lock, start, end, print_block, ingroup))

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
