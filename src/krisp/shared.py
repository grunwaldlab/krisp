import math
import tempfile
import itertools
from pathlib import Path
from .Amplicon import Amplicon, ConservedEndAmplicons


def prettyTime(time):
    """ Converts seconds into a readable minutes-seconds output format.

    Parameters
    ----------
    time : float
        Input time in seconds to convert

    Returns
    -------
    str
        A pretty version of the time to print
    """
    if time < 60:
        # Print seconds to 2 decimal places
        return f"{time:.2f} seconds"
    else:
        # Determine the number of minutes and remaining seconds
        minutes = int(time / 60)
        seconds = math.ceil(time - 60*minutes)
        # Determine if minutes and seconds are plural
        pl_min = "s" if (minutes > 1) else ""
        pl_sec = "s" if (seconds > 1) else ""
        return f"{minutes} minute{pl_min} and {seconds} second{pl_sec}"


def basename(filename):
    """Return the basename of a compressed fasta file.

    Parameters
    ----------
    filename : str
        The name of a file

    Returns
    -------
    str
        The basename of the file without fasta and compression endings

    """
    # Use Path library to get the filename from a full directory path and split
    filename = Path(filename).name.split('.')
    # Iterate and remove extensions from list
    omit = ["gz", "bz2", "fna", "fasta", "fa", "ffn", "frn"]
    while filename[-1] in omit:
        filename.pop()
    # Return what remains as the result
    return '.'.join(filename)


def simplename(filename):
    """Return the basename of a fasta file with all extensions removed.

    Parameters
    ----------
    filename : str
        The name of a file

    Returns
    -------
    str
        The basename of the file with no extensions

    """
    # Return what remains as the result
    return basename(filename).split('.')[0]


def tmpFile(workdir=None):
    """ Create a temporary file that won't be deleted on close.

    Parameters
    ----------
    workdir : str, optional
        Directory to store the newly created temporary file

    Returns
    -------
    str
        Filename of temporary file created
    """
    # Get name of working directory
    if workdir is None:
        workdir = tempfile.gettempdir()

    # Make a temporary file in a context s.t. it is closed automatically
    filename = ""
    with tempfile.NamedTemporaryFile(dir=workdir, delete=False) as fptr:
        # Get name parameter
        filename = fptr.name
    # Return filename
    return filename


def fileSize(fptr):
    """ Determines the length of an open file in bytes.

    This function determines the length of a file by moving the file pointer
    to the end of the file to get the position of the last byte. The pointer
    is returned to its original location unmodified.

    Parameters
    ----------
    fptr : open file-ptr like object

    Results
    -------
    int
        Size of the open file in bytes

    """
    # Get current pos
    origin = fptr.tell()

    # Move ptr to end of file
    fptr.seek(0, 2)

    # Get byte value as size
    file_size = fptr.tell()

    # Move ptr back to start and return size
    fptr.seek(origin)
    return file_size


def splitFilePtrs(filename, num_sections):
    """ Split a file into roughly equal sections.

    This function attempts to split a file into equally sized sections and
    return the start and end byte values for each section. The sections are
    only roughly equal as line lengths may differ and the file may not be
    equally divisable. The general algorithm is as follows:
        1) Open file for reading
        2) Get file size in bytes from open file pointer
        3) Determine the approximate size of each section in bytes
        4) Iterate to get official start and end bytes under the conditions:
            a) start byte is the beginning of a line


    Parameters
    ----------
    filename : str
        Name of file to get split positions

    num_sections : int
        The number of sections the file should be split into

    Returns
    -------
    list[tuple(int, int)]
        A list of start and end bytes which define each section

    """
    with open(filename, "r") as fptr:
        # Get size of file in bytes
        file_size = fileSize(fptr)

        # Divide size by num_sections to get aprox section sizes
        section_size = file_size // num_sections

        # Iterate through sizes and store pairs of start, end
        pairs = []
        start_byte = 0
        for i in range(num_sections):
            # Get approx end byte
            approx_end = start_byte + section_size
            if approx_end > file_size:
                approx_end = file_size

            # Move fileptr to approx end and then to start of line
            fptr.seek(approx_end)

            # Readline to get to next full line
            junk = fptr.readline()

            # Now iterate until the next line starts with a diff seq than this
            this_pos = fptr.tell()
            this_line = fptr.readline()
            next_pos = fptr.tell()
            next_line = fptr.readline()
            while (this_line.split(',')[0] == next_line.split(',')[0]):
                # iterate
                this_pos = next_pos
                this_line = next_line
                next_pos = fptr.tell()
                next_line = fptr.readline()
                # Check end of file
                if this_pos == next_pos:
                    break

            # Now this_line is different than next line, so next pos is end
            end_byte = next_pos

            # Save this value as end_byte
            pairs.append((start_byte, end_byte))

            # Print for testing
            fptr.seek(start_byte)
            start_byte = end_byte
        return pairs


def simplifyStream(stream):
    """ Simplify a stream by merging same sequences

    Parameters
    ----------
    stream : generator of Amplicon or ConservedEndAmplicons
        An iterable of Amplicons or ConservedEndAmplicons to simplify

    """
    # Iterate through stream and merge conserved sequences
    prev_item = None
    for item in stream:
        if item is None:
            # End of stream
            if prev_item is not None:
                yield prev_item
            return
        elif prev_item is None:
            # Just set as previous
            prev_item = item
        elif item == prev_item:
            # Same as previous so merge
            prev_item += item
        else:
            # Different so yield + reset
            yield prev_item
            prev_item = item

    # Yield last tuple if not None
    if prev_item is not None:
        yield prev_item


def removeSingletons(stream):
    """ Remove singletons from a stream

    Parameters
    ----------
    stream : generator of Amplicon or ConservedEndAmplicons
        An iterable of Amplicons or ConservedEndAmplicons to simplify

    Yields
        An iterable with all singletons removed

    """
    # Iterate through stream and merge conserved sequences
    prev_item = None
    singleton = True
    for item in stream:
        if item is None:
            # End of stream
            if prev_item is not None:
                yield prev_item
            return
        elif prev_item is None:
            # Just set as previous
            prev_item = item
            singleton = True
        elif item == prev_item:
            # Same as previous so yield
            yield prev_item
            prev_item = item
            singleton = False
        else:
            # Different so yield + reset
            if not singleton:
                yield prev_item
            prev_item = item
            singleton = True

    # Yield last item if not singleton
    if (prev_item is not None) and (not singleton):
        yield prev_item


def mergeSortedStreams(stream0, stream1):
    """ Merge two streams, maintaining sorted order

    Parameters
    ----------
    stream0 : Generator
        Generator of Amplicons or ConservedEndAmplicons

    stream1 : Generator
        Generator of Amplicons or ConservedEndAmplicons

    Yields
    ------
    Amplicons or ConservedEndAmplicons
        An iterable of sorted amplicons or ConservedEndAmplicons

    """
    # Pad streams with Nones so they don't end
    stream0 = itertools.chain(stream0, itertools.repeat(None))
    stream1 = itertools.chain(stream1, itertools.repeat(None))

    # Iterate while both are not None
    next_s0 = next(stream0)
    next_s1 = next(stream1)
    while (next_s0 is not None) or (next_s1 is not None):
        # Split on ordering
        if (next_s0 is None) or (next_s1 is not None) and (next_s1 < next_s0):
            # next_s1 is smaller
            yield next_s1
            next_s1 = next(stream1)
        else:
            # next_s0 is smaller
            yield next_s0
            next_s0 = next(stream0)


def intersectSortedStreams(stream0, stream1):
    """ Merge two ConservedEndAmplicons streams

    Parameters
    ----------
    stream0 : Generator
        Generator of ConservedEndAmplicons

    stream1 : Generator
        Generator of ConservedEndAmplicons

    Yields
    ------
    ConservedEndAmplicons
        An iterable of ConservedEndAmplicons

    """

    # Make sure each stream is simplified, i.e. duplicate elements are merged
    stream0 = simplifyStream(stream0)
    stream1 = simplifyStream(stream1)

    # Merge streams together
    merged = mergeSortedStreams(stream0, stream1)

    # Remove singletons and yield
    yield from simplifyStream(removeSingletons(merged))


def ampliconStream(filename, start=None, end=None):
    """ Converts a kmer file into a stream of Amplicons

    Parameters
    ----------
    filename : str
        Name of kmer file to read

    start : int, optional
        Starting byte of file to read

    end : int, optional
        End byte of file to read

    Yields
    ------
    Amplicon
        A stream of amplicons derived from the file with all duplicates merged

    """
    # Create a helper function to read from file
    def helper(filename, start, end):
        # Get tag from filename
        tag = simplename(filename)

        # Open file for reading
        with open(filename) as fptr:
            # Get start and end bytes
            if start is None:
                start = 0
            if end is None:
                end = fileSize(fptr)

            # Move file pointer to start
            fptr.seek(start)
            # Iterate until end is reached
            while fptr.tell() < end:
                # Get next line or break
                line = fptr.readline()
                if not line:
                    break
                # Convert line into an amplicon
                yield Amplicon.read(line, tag)

    # First get the raw amplicon stream
    stream = helper(filename, start, end)

    # Now simplify stream and yield
    yield from simplifyStream(stream)


def writeAmpliconStream(stream, output, mode="w"):
    """ Write an Amplicon stream to a file

    Parameters
    ----------
    stream : Iterable of amplicons
        An iterable of Amplicons to write

    output : str
        An output filename to write stream to

    """
    # Write all amplicons to output
    with open(output, mode) as fout:
        for amplicon in stream:
            amplicon.write(fout)


def writeAlignmentStream(stream, output, mode="w"):
    """ Write an alignment stream to a file

    Parameters
    ----------
    stream : Iterable of alignments
        An iterable of ConservedEndAmplicons to write

    output : str
        An output filename to write stream to

    """
    # Get amplicon stream
    def helper(align_stream):
        for align in align_stream:
            for ampl in align.amplicons:
                yield ampl

    ampl_stream = helper(stream)
    writeAmpliconStream(ampl_stream, output, mode=mode)



def alignmentStream(kmerfile, start=None, end=None, ingroup=None):
    """ Read a kmer file from start to end and yield alignments

    Parameters
    ----------
    filename : str
        Name of kmer file to read

    start : int, optional
        Starting byte of file to read

    end : int, optional
        End byte of file to read

    Yields
    ------
    Amplicon
        A stream of alignments derived from the file

    """
    # Read the kmerfile from start to end
    alignment = ConservedEndAmplicons(ingroup)
    for amplicon in ampliconStream(kmerfile, start, end):
        # Check if this amplicon can't fit into the alignment
        if not alignment.canAdd(amplicon):
            # Different so yield this alignment and reset
            yield alignment
            alignment = ConservedEndAmplicons(ingroup)
        # Add to alignment
        alignment += amplicon

    # Print last alignment if not empty
    if len(alignment):
        yield alignment


