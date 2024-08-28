import itertools
import subprocess
import multiprocessing
import shutil
import sys
import os
import time
from colorama import Fore, Back, Style
from .shared import *


def bisect(filename, key):
    """ Find the largest pointer value s.t. ptr.readline() < key

    Parameters
    ----------
    filename : str
        Name of file to bisect
    key : str
        Find the largest pointer value s.t. ptr.readline() < key

    Returns
    -------
    int
         Byte index for start of line

    """
    # Open filename for reading
    with open(filename, "r") as fptr:
        # Get size of file in bytes
        file_size = fileSize(fptr)

        # Check for Null input
        if key == "":
            return file_size, ""

        # bisect
        last_success = 0
        bracket = (0, file_size)
        prev_bracket = None
        while (bracket != prev_bracket):
            # Get starting values
            start, end = bracket

            # Get midpoint
            midpoint = (start + end) // 2

            # Move ptr to next full line
            fptr.seek(midpoint)
            junk = fptr.readline()
            midpoint = fptr.tell()

            # Update start, end
            line = fptr.readline()
            line = line.split(',')[0]
            if line < key:
                start = midpoint
                last_success = start
            else:
                end = midpoint

            # Update bracket
            prev_bracket = bracket
            bracket = (start, end)
        fptr.seek(last_success)
        return last_success, fptr.readline()


def mergeKmerFiles(filename0, filename1, output,
                   f0_start=None, f0_end=None, f1_start=None):
    """ Merge files in an output file

    Parameters
    ----------
    filename0 : str
        Name of first kmer file to merge
    filenam1 : str
        Name of second kmer file to merge
    output : str
        Name of output file to write
    f0_start : int, optional
        Start byte of first kmer file
    f0_end : int, optional
        End byte of first kmer file
    f1_start : int, optional
        Start byte of second kmer file

    Returns
    -------
    None
        Output file is written, None is returned
    """
    # generate alignment streams
    stream0 = alignmentStream(filename0, f0_start, f0_end)
    stream1 = alignmentStream(filename1, f1_start)

    # Take intersection of alignments
    alignments = intersectSortedStreams(stream0, stream1)

    # Write alignments back to file
    writeAlignmentStream(alignments, output)


def mergeKmerSingleCore(filename0, filename1, output):
    """ Merge two kmer files to output file in parallel

    Parameters
    ----------
    filename0 : str
        Name of first kmer file to merge
    filename1 : str
        Name of second kmer file to merge
    output : str
        Name of output file to write
    workdir : str, optional
        Name of work directory to write temporary files

    Returns
    -------
    None
        Output file is written, None is returned
    """
    # Now merge these files together using cat
    with open(output, "w") as fout:
        # Merge
        mergeKmerFiles(filename0, filename1, output)
        # Return merged filename
        return output


def mergeKmerParallel(filename0, filename1, output,
                      parallel, workdir=None):
    """ Merge two kmer files to output file in parallel

    Parameters
    ----------
    filename0 : str
        Name of first kmer file to merge
    filename1 : str
        Name of second kmer file to merge
    output : str
        Name of output file to write
    workdir : str, optional
        Name of work directory to write temporary files

    Returns
    -------
    None
        Output file is written, None is returned
    """
    # Start by making a list of temporary files
    tmp_files = [tmpFile(workdir) for i in range(parallel)]

    # Get the start and end iteration points for filename0
    f0_start_end = list(splitFilePtrs(filename0, parallel))

    # Find the bisection starting point for filename1
    f1_start = []
    with open(filename0) as fptr:
        for start, end in f0_start_end:
            # Get line at start
            fptr.seek(start)
            line = fptr.readline().split(',')[0]

            # Bisect on line and add pos to other_start
            pos, line2 = bisect(filename1, line)
            f1_start.append(pos)

    # Create a list of job args
    job_args = []
    for x, y, z in zip(tmp_files, f0_start_end, f1_start):
        job_args.append((filename0, filename1, x, *y, z))

    # Run jobs in parallel pool
    with multiprocessing.Pool(parallel) as pool:
        pool.starmap(mergeKmerFiles, job_args)

    # Now merge these files together using cat
    with open(output, "w") as fout:
        # Merge tmp_files into fptr.name
        for tmp in tmp_files:
            with open(tmp) as fin:
                fin.flush()
                shutil.copyfileobj(fin, fout)
                fout.flush()
        # Return merged filename
        return output


def mergeKmerJob(in_queue, out_queue, verbose=True):
    """ Process job for merging kmer files

    Parameters
    ----------
    in_queue : Multiprocessing Queue
        Input queue which holds jobs to run
    out_queue : Multiprocessing Queue
        Output queue to put finished jobs

    """
    # Get first job and iterate until None is encountered
    job_args = in_queue.get()
    while job_args is not None:
        # Split into args
        filename0, filename1, output, parallel, workdir = job_args

        # Print verbosity
        start_t = None
        if verbose:
            message = (f"Merging: {filename0} + {filename1} -> {output} using"
                       f" {parallel} cores")
            print(message, file=sys.stderr)
            start_t = time.time()

        # Merge files
        # Note that mergeKmerSingleCore is a temporary replacement for mergeKmerParallel which does not produce consistent results currently. See issue 17 on github for more details.
        #mergeKmerParallel(filename0, filename1, output, parallel, workdir)
        mergeKmerSingleCore(filename0, filename1, output)

        # Print verbosity
        if verbose:
            end_t = time.time()
            end_message = (f"=> Merged {filename0} + {filename1} -> {output}"
                           f" in {prettyTime(end_t-start_t)}")
            print(Fore.GREEN + end_message + Style.RESET_ALL, file=sys.stderr)

        # Add to output queue and get next job
        out_queue.put(output)
        job_args = in_queue.get()


def mergeFiles(files, output, parallel=1, workdir=None, verbose=True):
    """ Iteratively merges files until only one file is left

    Parameters
    ----------
    files : [str, str, ...]
        List of kmer filenames to merge
    output : str
        Ouptut file to write results
    parallel : int, optional
        Number of cores to use on this job
    wordir : str, optional
        Work directory to write temporary files
    verbose : bool, optional
        Print progress to stderr

    Returns
    -------
    None
        Output file is written
    """
    # Iterate through files and merge until only one file remains
    tmp_files = []
    counter = 0
    while len(files) > 1:
        # Get number of jobs to run and the number of cores per job
        jobs = len(files) // 2
        job_cores = []
        if parallel < jobs:
            job_cores = [1] * jobs
        else:
            for i in range(jobs):
                np = (parallel // jobs) + (i < parallel % jobs)
                job_cores.append(np)

        # Initialize processors with input and output queue
        processes = []
        job_queue = multiprocessing.Queue()
        fin_queue = multiprocessing.Queue()
        for i in range(jobs):
            p = multiprocessing.Process(target=mergeKmerJob,
                                        args=(job_queue, fin_queue, verbose))
            p.start()
            processes.append(p)

        # Get pairs of files to merge
        file_pairs = []
        while len(files) > 1:
            file_pairs.append((files.pop(), files.pop()))

        # Add jobs to queue
        for f0_f1, np in zip(file_pairs, job_cores):
            # Create temporary output file
            tmp_file = tmpFile(workdir)

            # Add job to queue
            f0, f1 = f0_f1
            job_queue.put((f0, f1, tmp_file, np, workdir))

        # now wait for jobs to finish and add to output list
        results = []
        for p in processes:
            result = fin_queue.get()
            results.append(result)
            tmp_files.append(result)

        # Send kill signal to processes and then join
        for p in processes:
            job_queue.put(None)

        # Call join on processes
        for p in processes:
            p.join()

        # Set new files to merge
        files = results + files

    # move result to output
    shutil.move(files[0], output)
