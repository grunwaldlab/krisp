from .shared import *


def filterAlignmentStream(stream, ingroup):
    """ Filter alignment stream based on ingroup

    Parameters
    ----------
    stream : generator of ConservedEndAmplicons
        A generator of alignments to filter
    ingroup : set
        A set-like object containing the simple names which denote the
        in-group, e.g. the set of genomes which should all contain a unique
        SNP which is identifiable

    Yields
    ------
    ConservedEndAlignments
        The filtered stream of alignments

    """
    # Iterate through every alignment in the stream
    for alignment in stream:
        # Set the ingroup
        alignment.setIngroup(ingroup)
        # Determine column numbers of unique SNP's : empty list -> None
        if len(alignment.ingroupUniqueColumns()) > 0:
            yield alignment


def filterAlignments(kmerfile, output, ingroup):
    # Read in alignments
    alignments = alignmentStream(kmerfile)

    # Filter alignment stream if ingroup is not empty
    if len(ingroup):
        alignments = filterAlignmentStream(alignments, ingroup)

    # Write to output
    writeAlignmentStream(alignments, output)
