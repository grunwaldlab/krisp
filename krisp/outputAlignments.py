import gzip
import multiprocessing
import shutil
import os
from .shared import *


def safePrint(alignments, counter, lock):
    """ Function to safely print alignments to the terminal in parallel """
    # Acquire lock
    with lock:
        # Iterate through alignments and print, increment counter
        for align in alignments:
            print(f"> Alignment_{counter.value}")
            print(align)
            counter.value += 1


def printAlignments(kmerfile, counter, lock,
                    start=None, end=None, print_block=10):
    """ Helper function to print alignments from kmer file """
    # Get alignment stream
    alignments = alignmentStream(kmerfile, start, end)
    alignments_to_print = []
    for alignment in alignments:
        # Add to alignments list
        alignments_to_print.append(alignment)
        # Print
        if len(alignments_to_print) >= print_block:
            # Print alignments
            safePrint(alignments_to_print, counter, lock)
            alignments_to_print = []
    if len(alignments_to_print):
        safePrint(alignments_to_print, counter, lock)


def printAlignmentsParallel(kmerfile, parallel, print_block=10):
    """ Function to print algignments from a file in parallel

    Parameters
    ----------

    kmerfile : str
        The input kmerfile to convert to an alignment file

    parallel : int
        The number of processors to use

    Returns
    -------
    int
        The number of alignments written

    """

    # Get the start and end iteration points for kmerfile
    fptr_start_end = list(splitFilePtrs(kmerfile, parallel))

    # Create a counter and lock
    counter = multiprocessing.Value('i', 0)
    lock = multiprocessing.Lock()

    # Create a list of job args
    job_args = []
    for start, end in fptr_start_end:
        job_args.append((kmerfile, counter, lock, start, end, print_block))

    # Start a group of processes
    processes = []
    for job in job_args:
        p = multiprocessing.Process(target=printAlignments,
                                    args=job)
        p.start()
        processes.append(p)

    # Call join on processes
    for p in processes:
        p.join()
    return counter.value


def writeAlignments(kmerfile, output, pid, start=None, end=None):
    """ Helper function to write alignments from kmer file """
    # Get alignment stream
    alignments = alignmentStream(kmerfile, start, end)

    # Open for writing and write output
    with gzip.open(output, 'wt') as fout:
        found = 0
        for i, alignment in enumerate(alignments):
            # Print allignment
            fout.write(f"> process_{pid}_{i}\n")
            fout.write(f"{str(alignment)}\n")
            found += 1
        fout.flush()
    return found


def writeAlignmentsParallel(kmerfile, output, parallel, workdir=None):
    """ Function to write alignments from a kmer file in parallel

    Parameters
    ----------

    kmerfile : str
        The input kmerfile to convert to an alignment file

    output : str
        The output alignment file

    parallel : int
        The number of processors to use

    workdir : str (default None)
        The working directory to place temporary files

    Returns
    -------
    int
        The number of alignments written

    """

    # Get the start and end iteration points for kmerfile
    fptr_start_end = list(splitFilePtrs(kmerfile, parallel))

    # Create a list of temporary output files
    tmp_files = [tmpFile(workdir) for i in range(parallel)]

    # Create a list of job args
    job_args = []
    for i in range(parallel):
        start, end = fptr_start_end[i]
        pid = i
        tmp = tmp_files[i]
        job_args.append((kmerfile, tmp, pid, start, end))

    # Run processes in pool
    with multiprocessing.Pool(parallel) as pool:
        counts = pool.starmap(writeAlignments, job_args)

    # Merge tmpfiles into output
    if output.endswith('.gz'):
        with open(output, "wb") as fout:
            for tmp in tmp_files:
                with open(tmp, "rb") as fin:
                    shutil.copyfileobj(fin, fout)
    else:
        with open(output, "w") as fout:
            for tmp in tmp_files:
                with gzip.open(tmp, "rt") as fin:
                    shutil.copyfileobj(fin, fout)
    return sum(counts)
