#! /bin/env python3
import sys
import argparse
import fileinput
import itertools
import multiprocessing
import subprocess


# Dictionary of Watson-Crick complements for DNA
COMP_MAP = {'A': 'T', 'T': 'A', 'a': 't', 't': 'a',
            'G': 'C', 'C': 'G', 'g': 'c', 'c': 'g',
            'R': 'Y', 'Y': 'R', 'r': 'y', 'y': 'r',
            'M': 'K', 'K': 'M', 'm': 'k', 'k': 'm',
            'S': 'S', 'W': 'W', 's': 's', 'w': 'w',
            'B': 'V', 'V': 'B', 'b': 'v', 'v': 'b',
            'D': 'H', 'H': 'D', 'd': 'h', 'h': 'd',
            'N': 'N', 'n': 'n'}

# IUPAC nucleobase generalizations
IUPAC_BASE = {'R': ['A', 'G'],
              'Y': ['C', 'T'],
              'S': ['G', 'C'],
              'W': ['A', 'T'],
              'K': ['G', 'T'],
              'M': ['A', 'C'],
              'B': ['C', 'G', 'T'],
              'D': ['A', 'G', 'T'],
              'H': ['A', 'C', 'T'],
              'V': ['A', 'C', 'G'],
              'N': ['A', 'C', 'G', 'T'],
              'r': ['a', 'g'],
              'y': ['c', 't'],
              's': ['g', 'c'],
              'w': ['a', 't'],
              'k': ['g', 't'],
              'm': ['a', 'c'],
              'b': ['c', 'g', 't'],
              'd': ['a', 'g', 't'],
              'h': ['a', 'c', 't'],
              'v': ['a', 'c', 'g'],
              'n': ['a', 'c', 'g', 't']}


def sortPipe(np=None, mem=None, cols=None):
    """ Create a sorting subprocess Pipe

    Parameters
    ----------
    np : int (optional)
        The number of cores to use for sorting

    mem : str (optional)
        The amount of memory to use during sorting. See linux sort '-S'

    cols : list of int (optional)
        Columns to sort on, assumes 0-indexing and ',' seperating

    Returns
    --------
    subprocess Pipe
        An active subprocess pipe that can be written to and read from

    """
    # Get default values
    command = f"LC_ALL=C sort"
    if mem is not None:
        command += f" -S {mem}"
    if np is not None:
        command += f" --parallel={np}"
    if (cols is not None):
        command += f" -t,"
        for c in cols:
            command += f" -k{c+1},{c+1}"
    # Start subprocess and return
    return subprocess.Popen(command,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            text=True,
                            shell=True)


def sortInPlace(filename, np=None, mem=None, cols=None):
    """ Sorts a file in place

    Parameters
    ----------
    filename : str
        The file to sort

    np : int (optional)
        The number of cores to use for sorting

    mem : str (optional)
        The amount of memory to use during sorting. See linux sort '-S'

    cols : list of int (optional)
        Columns to sort on, assumes 0-indexing and ',' seperating

    Returns
    --------
    None
        The input file is sorted in place

    """

    # Get default values
    command = f"LC_ALL=C sort {filename} -o {filename}"
    if mem is not None:
        command += f" -S {mem}"
    if np is not None:
        command += f" --parallel={np}"
    if (cols is not None):
        command += f" -t,"
        for c in cols:
            command += f" -k{c+1},{c+1}"
    # Run subprocess
    process = subprocess.Popen(command, shell=True)
    process.wait()


class kstream:
    """
    A highly flexible class to read and parse kmers from a fasta file. The
    class operates on the idea of filtering and manipulating kmer streams. The
    input fasta file is converted to a kmer stream which is then manipulated
    based on the passed arguments.
    """

    def __init__(self, sequences=None, kmers=None, complements=False,
                 canonicals=False, allow=None, disallow=None, omitsoft=False,
                 mapsoft=False, expandiupac=False, split=None, sort=False,
                 sortmem=None, sortcols=None, sortnp=1, parallel=1):
        """
        Parameters
        ----------
        sequences : iterable of str, or fasta filename
            An iterable of sequence strings to parse or a fasta file

        kmers : int or [int]
            The desired kmer length or a list of desired kmer lengths

        complements : bool (default False)
            True will add the reveres complements to the resulting stream

        canonicals : bool (default False)
            True will convert all kmers to canonical form (alphabetical order)

        allow : str
            A list of allowed characters in kmers. Kmers will be omitted if
            they contain any character outside of allowed, e.g allow='ATGC'
            will only allow kmer with an 'A', 'T', 'G', 'C'.

        disallow : str
            A list of disallowed characters in kmers. The opposite of allow.

        omitsoft : bool (default False)
            Disallow any kmer containing a lowercase character.

        mapsoft : bool (default False)
            Map lowercase characters to uppercase in all kmers.

        expandiupac : bool (default False)
            Expand any IUPAC characters and yield each as a separate kmer.

        split : [int]
            A list of positions to split the resulting kmers on to create a
            series of ',' seperated sub-kmers.

        sort : bool (default False)
            Sort the resulting character stream. Requires linux sort utility.
            See sortmem, sortcols, sortnp for additional options.

        sortmem : str
            The amount of memory to delegate to the sorting routine. See sort
            linux command. Can be amount: "1G", or a %: "50%".

        sortcols : [int]
            A list of integer positions (0-indexed) to sort on. Default sorts
            on entire kmer. Useful with split command to sort on 1st, 3rd
            columns, e.g. sortcols=[0,2].

        sortnp : int
            Number of cores to dedicate to sorting

        parallel : int
            Number of cores to dedicate to kmer extraction

        Returns
        -------
        kstream instance
            Any kstream instance can be iterated through as though it was a
            list, assuming sequences was not None when initiated. Or, a kstream
            instance can be initiated and called repeatedly as a functor to
            retrieve an iterable instance.

        """

        # Parse stream options and add to parsers
        # Each parser is a function which takes an input stream, modifies it,
        # And yields the result. For example, self._complements will take a
        # Stream of kmers, and yield every input kmer, and its complement.
        self.parsers = []
        if kmers is not None:
            if isinstance(kmers, int):
                self.kmers = [kmers]
            else:
                self.kmers = list(kmers)
            self.parsers.append(self._kmers)
        if omitsoft is True:
            if mapsoft is True:
                raise ValueError("can't omit and map soft masked nucleotides")
            self.parsers.append(self._omitsoft)
        if mapsoft is True:
            self.parsers.append(self._mapsoft)
        if complements is True:
            if canonicals is True:
                raise ValueError("canonicals conflicts with complements")
            self.parsers.append(self._complements)
        if allow is not None:
            self.allow = set(allow)
            self.parsers.append(self._allow)
        if disallow is not None:
            self.disallow = set(disallow)
            self.parsers.append(self._disallow)
        if expandiupac is True:
            self.parsers.append(self._expandiupac)
        if canonicals is True:
            self.parsers.append(self._canonicals)
        if split is not None:
            if isinstance(split, int):
                self.split = [split]
            else:
                self.split = list(split)
            self.parsers.append(self._split)

        # Check if sorting is specified
        self.sort = sort
        if self.sort is True:
            self.sortnp = sortnp
            self.sortmem = sortmem
            self.sortcols = sortcols

        # Set parallel flag
        self.parallel = parallel

        # Store sequences for future call to iter
        self.sequences = sequences

    def write(self, filename, sequences=None):
        """ Write sequences to a file

        Parameters
        ----------
        filename : str
            Name of file to write kmers to

        sequences : [str] or str (Default None)
            The sequences to write, if None, assumes kstream was instanciated
            with sequences.

        Returns
        -------
        int
            The number of kmers written to the file

        """
        # Parse sequences
        if sequences is None:
            sequences = self.sequences
        sequences = self._parsed_kmers(sequences)

        # Check if RNA
        is_rna, sequences = self._detect_RNA(sequences)

        # Map to DNA if input is RNA
        if is_rna:
            sequences = self._to_DNA(sequences)

        # Keep count of sequences found
        count = 0

        # Open file for writing
        with open(filename, 'w') as fout:
            # Split on parallel
            if self.parallel == 1:
                # Apply filters
                for f in self.parsers:
                    sequences = f(sequences)

                # Map back to RNA
                if is_rna:
                    sequences = self._to_RNA(sequences)

                # Write sequences
                for seq in sequences:
                    print(f"{seq}", file=fout)
                    count += 1
            else:
                # Run parsing in parallel using a pool of workers
                with multiprocessing.Pool(self.parallel) as pool:
                    # Parse sequences using pools imap, result is nested list
                    for result in pool.imap_unordered(self._parallel_job,
                                                      sequences):
                        # Map back to RNA
                        if is_rna:
                            result = self._to_RNA(result)

                        # Write to file
                        for seq in result:
                            print(f"{seq}", file=fout)
                            count += 1
            # Make sure that file contents have been written
            fout.flush()

        # Check if sort
        if self.sort:
            # Sort file in place
            sortInPlace(filename,
                        np=self.sortnp,
                        mem=self.sortmem,
                        cols=self.sortcols)

        # Return number of kmers found
        return count

    def __call__(self, sequences):
        """ Convert sequences into a kmer stream

        Parameters
        ----------

        sequences : iterable
            The iterable of sequences or a fasta filename

        Yields
        ------
        str
            The resulting kmers

        """
        # Parse sequences
        sequences = self._parsed_kmers(sequences)

        # Check if RNA
        is_rna, sequences = self._detect_RNA(sequences)

        # Map to DNA if input is RNA
        if is_rna:
            sequences = self._to_DNA(sequences)

        # Split on parallel
        if self.parallel == 1:
            # Apply filters
            for f in self.parsers:
                sequences = f(sequences)

            # Check if sort
            if self.sort:
                # Create sort pipe
                sort_pipe = sortPipe(np=self.sortnp,
                                     mem=self.sortmem,
                                     cols=self.sortcols)
                # Write sequences to sort
                for seq in sequences:
                    print(f"{seq}", file=sort_pipe.stdin)
                sort_pipe.stdin.close()

                # set values back into sequences
                sequences = (line.strip() for line in sort_pipe.stdout)

            # Map back to RNA
            if is_rna:
                sequences = self._to_RNA(sequences)

            # Yield sequences
            yield from sequences
        else:
            # Run parsing in parallel using a pool of workers
            sort_pipe = None
            if self.sort:
                sort_pipe = sortPipe(np=self.sortnp,
                                     mem=self.sortmem,
                                     cols=self.sortcols)
            with multiprocessing.Pool(self.parallel) as pool:
                # Parse sequences using pools imap, results is a nested list
                for result in pool.imap_unordered(self._parallel_job,
                                                  sequences):

                    # Map back to RNA
                    if is_rna:
                        result = self._to_RNA(result)

                    # Check if sort
                    if self.sort:
                        # Write sequences to sort
                        for seq in result:
                            print(f"{seq}", file=sort_pipe.stdin)
                    else:
                        yield from result
            # Yield now if sort
            if self.sort:
                # Close pipe and yield result
                sort_pipe.stdin.close()
                yield from (line.strip() for line in sort_pipe.stdout)

    def _parallel_job(self, seq):
        """ Convert a single seq and return a list

        Parameters
        ----------
        seq : str
            A single input string to convert into kmers

        Returns
        -------
        [str]
            A list of kmer strings

        """
        sequences = (seq,)
        for f in self.parsers:
            sequences = f(sequences)
        return list(sequences)

    def __iter__(self):
        """ Return iterator of kmers from sequences provided at creation """
        return iter(self.__call__(self.sequences))

    def _parsed_kmers(self, sequences):
        """ Utility function to extract sequences from a file / iterable

        Parameters
        ----------
        sequences : [str] or str
            The iterable of sequences or a fasta filename

        Yields
        ------
        str
            An iterable of kmers

        """
        # Check if sequences is a file or an iterable
        if isinstance(sequences, str):
            # Filename
            sequences = self._read_file(sequences)

        # Check if fasta format and parse stream
        is_fasta, sequence = self._detect_FASTA(sequences)
        if is_fasta:
            sequences = self._parse_FASTA(sequences)
        else:
            sequences = self._parse_seqs(sequences)

        yield from sequences

    def _read_file(self, filename):
        """ Return iterator to file lines

        Parameters
        ----------
        filename : str
            The input Fasta file to read

        Yields
        ------
        str
            The lines of the input file

        """
        for line in fileinput.input(filename,
                                    openhook=fileinput.hook_compressed):
            # Try to decode
            try:
                tmp = line.decode()
                yield tmp
            except (UnicodeDecodeError, AttributeError):
                yield line

    def _detect_RNA(self, sequences):
        """ Checks if sequences is RNA or DNA, sets value and re-yields

        Parameters
        ----------
        sequences : Iterable of str
            The iterable of sequences

        Returns
        -------
        tuple of bool, iterable of str
            bool is True if sequences are RNA. Second return value is the
            Unmodified iterable of input sequences.

        """
        found = []
        is_rna = None
        for seq in sequences:
            # Add seq to queue
            found.append(seq)
            if 'T' in seq or 't' in seq:
                is_rna = False
                break
            if 'U' in seq or 'u' in seq:
                is_rna = True
                break
        # return tuple
        return (is_rna, itertools.chain(found, sequences))

    def _detect_FASTA(self, sequences):
        """ Checks if sequences is RNA or DNA, sets value and re-yields

        Parameters
        ----------
        sequences : generator of str
            The stream of lines to parse

        Returns
        ------
        tuple of bool, iterable of str
            bool is True if the sequence is from a Fasta file, else False
            The iterable is the unmodified sequence stream

        """
        lines = []
        is_fasta = None
        for seq in sequences:
            # Add seq to lines
            lines.append(seq)
            if '>' in seq:
                is_fasta = True
                break
            else:
                is_fasta = False
                break
        # return tuple
        return (is_fasta, itertools.chain(lines, sequences))

    def _parse_seqs(self, lines):
        """ Parses lines to remove spaces and newlines

        Parameters
        ----------
        lines : generator of str
            The stream of lines to parse

        Yields
        ------
        str
            The parsed stream of lines

        """
        for line in lines:
            yield line.strip()

    def _parse_FASTA(self, lines):
        """ Parse fasta file and yield sequences

        Parameters
        ----------
        lines : generator of str
            The stream of lines to convert into sequences

        Yields
        ------
        str
            A single sequence for each fasta entry

        """

        sequence = ""
        for line in lines:
            # Remove any newlines and spaces
            line = line.strip()
            if line.startswith('>'):
                if len(sequence):
                    yield sequence
                sequence = ""
            else:
                # Add line to sequence
                sequence += line
        if len(sequence):
            yield sequence

    def _to_DNA(self, sequences):
        """ Converts sequences to DNA

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        return (seq.replace('U', 'T').replace('u', 't') for seq in sequences)

    def _to_RNA(self, sequences):
        """ Converts sequences to RNA

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        return (seq.replace('T', 'U').replace('t', 'u') for seq in sequences)

    def _kmers(self, sequences):
        """ Convert sequecnes to kmers

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Dummy funtion
        def _extract(seq):
            # Iterate through every length
            for k in self.kmers:
                # Iterate through each index and yield
                for i in range(len(seq) - k + 1):
                    tmp = [seq[i + j] for j in range(k)]
                    yield ''.join(tmp)

        # Extract every kmer from every sequence and yield
        for seq in sequences:
            yield from _extract(seq)

    def _get_complement(self, seq):
        """ Returns complement to a single sequence

        Parameters
        ----------
        seq : str
            sequence to get reverse complement of

        Returns
        -------
        str
            The reverse complement

        """
        new_seq = ''.join([COMP_MAP[b] for b in reversed(seq)])
        return new_seq

    def _complements(self, sequences):
        """ Add complements to stream

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        for seq in sequences:
            yield seq
            yield self._get_complement(seq)

    def _canonicals(self, sequences):
        """ Only add canonical kmer to stream

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        for seq in sequences:
            yield min(seq, self._get_complement(seq))

    def _allow(self, sequences):
        """ Return allowed sequences

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Dummy function for filter
        def _isallowed(seq):
            return set(seq).issubset(self.allow)
        return filter(_isallowed, sequences)

    def _disallow(self, sequences):
        """ Omit disallowed sequences

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Dummy function for filter
        def _isallowed(seq):
            return set(seq).isdisjoint(self.disallow)
        return filter(_isallowed, sequences)

    def _omitsoft(self, sequences):
        """ Omit sequences containing softmasked nucleotides

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences


        """
        return filter(lambda x: x.isupper(), sequences)

    def _mapsoft(self, sequences):
        """ Map sequences containing softmasking to uppercase

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Map soft masked sequences to uppercase
        return (seq.upper() for seq in sequences)

    def _expandiupac(self, sequences):
        """ Expand iupac characters in sequences

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to parse

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Iterate through each sequence and call _expand
        for seq in sequences:
            # Find all IUPAC bases
            iupac_bases = []
            iupac_options = []
            for i, b in enumerate(seq):
                if b in IUPAC_BASE:
                    iupac_bases.append(i)
                    iupac_options.append(IUPAC_BASE[b])
            # Check IUPAC bases found:
            if len(iupac_bases):
                # Iterate through every combination of options
                seq = list(seq)
                for combo in itertools.product(*iupac_options):
                    # Place bases
                    for i, b in zip(iupac_bases, combo):
                        seq[i] = b
                    # yield
                    yield ''.join(seq)
            else:
                # No iupac so just yield
                yield seq

    def _split(self, sequences):
        """ Split sequence into ',' seperated columns

        Parameters
        ----------
        sequences : generator of str
            The stream of sequences to split

        Yields
        ------
        str
            The parsed stream of sequences

        """
        # Iterate through every sequence
        for seq in sequences:
            # Split seq into chunks
            pos_parts = []
            neg_parts = []
            for size in self.split:
                if size >= 0:
                    pos_parts.append(seq[:size])
                    seq = seq[size:]
                else:
                    neg_parts.append(seq[size:])
                    seq = seq[:size]
            # Merge parts and yield
            yield ','.join(pos_parts + [seq] + neg_parts)


def parseArgs(sys_args):
    """ Parse command line args """
    parser = argparse.ArgumentParser(
            description=("Read and parse kmers from fasta or kmer stream\n"
                         "Compatible with gz, bz2, and stdin."),
            prog="kstream",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
            "file",
            nargs="?",
            type=str,
            help="Fasta file to read. .gz, .bz2, default stdin",
            default="-")
    parser.add_argument(
            "-k",
            "--kmers",
            type=int,
            nargs="+",
            help="Convert sequences into kmers of given length(s).")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
            "--canonicals",
            action="store_true",
            help="Print canonical sequences (alphabetically first)")
    group.add_argument(
            "--complements",
            action="store_true",
            help="Add reverse complement to stream")
    parser.add_argument(
            "--disallow",
            type=str,
            help="Omit sequences containing dissallowed nucleotides")
    parser.add_argument(
            "--allow",
            type=str,
            help="Only accept sequences containing allowed nucleotides")
    parser.add_argument(
            "--expand-iupac",
            action="store_true",
            help="Expand IUPAC nucleotide codes (including N's)")
    parser.add_argument(
            "--omit-softmask",
            action="store_true",
            help="Omit sequences containing soft masking")
    parser.add_argument(
            "--map-softmask",
            action="store_true",
            help="Unmask sequences containing soft masking")
    parser.add_argument(
            "--split",
            nargs="+",
            type=int,
            help="Split kmers into columns and delimit by ','")
    parser.add_argument(
            "-p",
            "--parallel",
            type=int,
            default=1,
            help="Number of processors to use. Default 1")
    parser.add_argument(
            "-s",
            "--sort",
            action="store_true",
            help="Sort resulting kmers")
    parser.add_argument(
            "--sort-np",
            type=int,
            default=1,
            help="Number of processores to use for sorting")
    parser.add_argument(
            "--sort-mem",
            type=str,
            help="Amount of memory to use, see linux sort mem usage")
    parser.add_argument(
            "--sort-cols",
            nargs="+",
            type=int,
            help="Sort based on these columns, 0-based indexing")
    parser.add_argument(
            "--output",
            help="Write output to file as opposed to terminal")
    parser.add_argument(
            '--version',
            action='version',
            version='%(prog)s 1.0')

    # Parse args and return
    return parser.parse_args(sys_args)


def main():
    # Get command line args
    args = parseArgs(sys.argv[1:])

    # Create kstream instance
    streamer = kstream(kmers=args.kmers,
                       complements=args.complements,
                       canonicals=args.canonicals,
                       allow=args.allow,
                       disallow=args.disallow,
                       omitsoft=args.omit_softmask,
                       mapsoft=args.map_softmask,
                       expandiupac=args.expand_iupac,
                       split=args.split,
                       parallel=args.parallel,
                       sort=args.sort,
                       sortnp=args.sort_np,
                       sortmem=args.sort_mem,
                       sortcols=args.sort_cols)

    # Pass args.file into kstream and print or write result
    if args.output is not None:
        with open(args.output, "w") as fout:
            for seq in streamer(args.file):
                print(seq, file=fout)
    else:
        for seq in streamer(args.file):
            print(seq)


if __name__ == "__main__":
    main()
