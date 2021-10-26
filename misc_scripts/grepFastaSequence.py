from Bio import SeqIO
import sys
import gzip


def findAllSubstrings(string, substring):
    # Find first matching
    pos = string.find(substring)
    while (pos != -1):
        # Yield this finding
        yield pos
        # Get next position
        pos = string.find(substring, pos+1)


def findSequenceInFasta(fasta, seq, padding):
    # Iterate through every fasta entry in file
    fasta_entries = []
    if fasta.endswith(".gz"):
        fasta_entries = SeqIO.parse(gzip.open(fasta, "rt"), 'fasta')
    else:
        fasta_entries = SeqIO.parse(open(fasta, "r"), 'fasta')
    
    for fasta_entry in fasta_entries:
        # Get name and sequence
        name, sequence = fasta_entry.id, str(fasta_entry.seq)
        # Iterate through all positions containing seq
        for pos in findAllSubstrings(sequence, seq):
            start = int(max(0, pos - padding))
            end = int(pos + len(seq) + padding)
            yield (name, start + 1, end + 1, sequence[start: end])


def findComplementInFasta(fasta, seq, padding):
    # Get complement of seq
    mapping = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    rev_comp = ''.join(mapping[b] for b in reversed(seq))

    # Call find
    for name, start, end, sequence in findSequenceInFasta(fasta, rev_comp, padding):
        # Return sequence in lowercase to denote minus strand
        sequence = ''.join(mapping[b] for b in reversed(sequence))
        yield (name, start, end, sequence.lower())


def findInFasta(fasta, seq, padding):
    # Look for forward sequence, then reverse complement and yield
    for result in findSequenceInFasta(fasta, seq, padding):
        yield result
    for result in findComplementInFasta(fasta, seq, padding):
        yield result


if __name__ == "__main__":
    # Get fasta file, seq, and padding
    fasta = sys.argv[1]
    seq = sys.argv[2]
    padding = int(sys.argv[3])
    
    # Search for sequence in fasta and print all finding
    for name, start, end, sequence in findInFasta(fasta, seq, padding):
        print(f"{name.ljust(8)} {str(start).ljust(8)}   {sequence}   {end}")
