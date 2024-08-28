import sys
import re
import primer3
from statistics import mean
from Bio.Data import IUPACData
from prettytable import PrettyTable
# from pudb.remote import set_trace

UNKNOWN_CHAR = "?"
iupac_key = {tuple((x for x in sorted(v))): k for k, v in
             IUPACData.ambiguous_dna_values.items()}
iupac_key[(UNKNOWN_CHAR,)] = 'N'

primer3_col_names = [
    'PRIMER_PAIR_0_PRODUCT_SIZE',
    'PRIMER_PAIR_0_PENALTY',
    'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE',
    'PRIMER_LEFT_0_PENALTY', 'PRIMER_RIGHT_0_PENALTY',
    'PRIMER_LEFT_0_TM', 'PRIMER_RIGHT_0_TM',
    'PRIMER_LEFT_0_GC_PERCENT', 'PRIMER_RIGHT_0_GC_PERCENT',
    'PRIMER_LEFT_0_SELF_ANY_TH', 'PRIMER_RIGHT_0_SELF_ANY_TH',
    'PRIMER_LEFT_0_SELF_END_TH', 'PRIMER_RIGHT_0_SELF_END_TH',
    'PRIMER_LEFT_0_HAIRPIN_TH', 'PRIMER_RIGHT_0_HAIRPIN_TH',
    'PRIMER_LEFT_0_END_STABILITY', 'PRIMER_RIGHT_0_END_STABILITY',
    'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH',
]
primer3_col_key = {n: n.replace("PRIMER_", "").replace("_0", "").lower() for n in primer3_col_names}

# class ForkedPdb(pdb.Pdb):
#     """A Pdb subclass that may be used
#     from a forked multiprocessing child
# 
#     """
#     def interaction(self, *args, **kwargs):
#         _stdin = sys.stdin
#         try:
#             sys.stdin = open('/dev/stdin')
#             pdb.Pdb.interaction(self, *args, **kwargs)
#         finally:
#             sys.stdin = _stdin

def collapse_to_iupac(seqs):
    """Combine sequences into a consensus using IUPAC ambiguity codes

    Parameters:
    -----------
    seqs : list of str
        The sequences to combine.

    Returns
    -------
    str
        The combined sequence
    """
    seq_lens = [len(x) for x in seqs]
    max_len = max(seq_lens)
    if len(set(seq_lens)) != 1:  # TODO: replace with alignment
        return '-' * max_len
    output = []
    for i in range(max_len):
        column = {s[i] for s in seqs}
        if "*" in column or "N" in column or UNKNOWN_CHAR in column:
            output.append('N')
        else:
            output.append(iupac_key[tuple(sorted(column))])
    return "".join(output)


def parse_primer3_settings(file_path):
    """
    Reads primer3 BoulderIO format, assuming only global settings (starts with PRIMER_),
    and returns a dict in the format used by primer3-py.
    """
    def to_number_if_can(x):
        try:
            if int(float(x)) == float(x) and '.' not in x:
                return int(x)
            else:
                return float(x)
        except ValueError:
            return x

    with open(file_path) as handle:
        options = dict([tuple(l.strip().split('=')) for l in handle.readlines()])
        for opt, val in options.items():
            if ' ' in val or ';' in val:
                val = re.split('[ ;]+', val)
                val = [to_number_if_can(v) for v in val]
                if ',' in val or '-' in val[0]:
                    val = [[to_number_if_can(x) for x in re.split('[,-]+', v)] for v in val]
            elif ',' in val or '-' in val:
                val = re.split('[,\-]+', val)
                val = [to_number_if_can(v) for v in val]
            else:
                val = to_number_if_can(val)
            options[opt] = val
    return options

def _format_p3_output(p3_out):
    """Reformat data for best primer pair for CSV output"""
    return {primer3_col_key[n]: p3_out[n] for n in primer3_col_names}

def run_primer3(template, target_start, target_len,
                options=None,
                tm=(53, 68),
                gc=(40, 70),
                amp_size=(80, 300),
                primer_size=(25, 35),
                max_sec_tm=40,
                gc_clamp=1,
                max_end_gc=4):
    if options is None:
        global_options = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,  # 1 == True
            # 'PRIMER_PICK_INTERNAL_OLIGO': 1,  # 1 == True
            'PRIMER_PICK_RIGHT_PRIMER': 1,  # 1 == True
            'PRIMER_LIBERAL_BASE': 1,  # 1 == True
            'PRIMER_OPT_SIZE': mean(primer_size),
            'PRIMER_MIN_SIZE': primer_size[0],
            'PRIMER_MAX_SIZE': primer_size[1],
            # 'PRIMER_INTERNAL_MAX_SIZE': len(crrna_seq),
            'PRIMER_OPT_TM': mean(tm),
            'PRIMER_MIN_TM': tm[0],
            'PRIMER_MAX_TM': tm[1],
            'PRIMER_MIN_GC': gc[0],
            'PRIMER_MAX_GC': gc[1],
            'PRIMER_MAX_POLY_X': 4,  # The maximum allowable length of a mononucleotide repeat
            'PRIMER_MAX_NS_ACCEPTED': 0,  # The maximum number of Ns
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,  # 1 == True
            'PRIMER_MAX_SELF_ANY_TH': max_sec_tm,  # describes the tendency of a primer to bind to itself
            'PRIMER_MAX_SELF_END_TH': max_sec_tm,  # 3' binding to itself (primer dimers)
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': max_sec_tm,  # the tendency of the left primer to bind to the right primer
            'PRIMER_PAIR_MAX_COMPL_END_TH': max_sec_tm,  # primer heterodimers
            'PRIMER_MAX_HAIRPIN_TH': max_sec_tm,  # hairpins
            'PRIMER_PRODUCT_SIZE_RANGE': [amp_size],
            'PRIMER_GC_CLAMP': gc_clamp,
            'PRIMER_MAX_END_GC': max_end_gc,
        }
    else:
        global_options = parse_primer3_settings(options)

    p3_output = primer3.bindings.design_primers(
        {
            'SEQUENCE_TEMPLATE': "".join(template),
            # 'SEQUENCE_INTERNAL_OLIGO': "".join(crrna_seq),
            'SEQUENCE_TARGET': [target_start, target_len]
        },
        global_options
    )
    return p3_output


class Amplicon:
    def __init__(self, primer, diagnostic, reverse, *labels):
        """ Store values """
        # Store sequences
        self.primer = primer
        self.diagnostic = diagnostic
        self.reverse = reverse

        # Store labels as a list of filenames
        self.labels = sorted(labels)

    @property
    def sequence(self):
        """ Return the sequence for this amplicon """
        return f"{self.primer}{self.diagnostic}{self.reverse}"

    def _labelsToString(self):
        """ Helper function to convert labels into a string """
        # Convert labels to a dictionary for counting
        counts = {}
        for label in self.labels:
            if label not in counts:
                counts[label] = 1
            else:
                counts[label] += 1

        # Merge all labels into a single label list
        label_list = []
        for filename, count in sorted(counts.items()):
            if count == 1:
                label_list.append(f"{filename}")
            else:
                label_list.append(f"{filename}({count})")
        return ';'.join(label_list)

    def _stringToLabels(string):
        """ Helper function to convert a string into a series of labels """
        # Convert string to a list of labels
        labels = []

        # Read through tags to get label
        for label in string.split(';'):
            # Strip whitespace
            label = label.strip()

            # Check for multiplier
            if '(' in label:
                name, multiplier = label.split('(')
                multiplier = int(multiplier.strip(')'))
                labels += [name] * multiplier
            else:
                labels += [label]
        return labels

    def __str__(self):
        """ Returns a string representation of this amplicon """
        return f"{self.sequence} : {self._labelsToString()}"

    def __eq__(self, other):
        """ Returns true if same sequence as other

        Parameters
        ----------
        other : Amplicon
            The other amplicon to compare with

        Returns
        -------
        bool
            True if self and other have the same sequence, else False

        Raises
        ------
        TypeError
            If other is not an Amplicon instance

        """
        # Check Valid type
        if not isinstance(other, Amplicon):
            raise TypeError(f"Can't compare Amplicon with type {type(other)}")
        return self.sequence == other.sequence

    def __add__(self, other):
        """ Merges two amplicons and their labels

        Parameters
        ----------
        other : Amplicon
            The other amplicon to add

        Returns
        -------
        Amplicon
            A new Amplicon instance containing the merged labels

        Raises
        ------
        TypeError
            If other is not an Amplicon instance
        ValueError
            If the amplicons have different sequences

        """
        # Check Valid type
        if not isinstance(other, Amplicon):
            raise TypeError(f"Can't compare Amplicon with type {type(other)}")

        # Check valid addition
        if not self.__eq__(other):
            raise ValueError("Can't merge Amplicons with different sequences")

        # Create a new amplicon to return
        return Amplicon(self.primer,
                        self.diagnostic,
                        self.reverse,
                        *self.labels,
                        *other.labels)

    def __lt__(self, other):
        """ Returns true if self comes before other in terms of primer ordering

        Parameters
        ----------
        other : Amplicon
            The other amplicon to compare with

        Returns
        -------
        bool
            True if self comes before other in terms of primer sequences

        Raises
        ------
        TypeError
            If other is not an Amplicon instance

        """
        # Check Valid type
        if not isinstance(other, Amplicon):
            raise TypeError(f"Can't compare Amplicon with type {type(other)}")

        # Check primers
        return (self.primer, self.reverse) < (other.primer, other.reverse)

    def read(string, filename):
        """ Class method to read a string from a file to an Amplicon

        Parameters
        ----------
        string : str
            An input string to create an amplicon. Assumes the format:
            primer, diagnostic, reverse, tag0; tag1; tag2

        filename : str
            The default label to add to the Amplicon

        """
        # Split string by commas
        fields = string.strip().split(',')
        if len(fields) in [3, 4]:
            # Split fields
            primer, diag, reverse, *tags = fields

            # Format tags
            labels = []
            if len(tags) == 0:
                # Set label to be filename
                labels = [filename]
            else:
                labels = Amplicon._stringToLabels(tags[0])
            # Create amplicon and return
            return Amplicon(primer, diag, reverse, *labels)
        else:
            raise ValueError(f"Unrecognised string format : {string}")
            #removed "return None" as it was unreachable

    def write(self, fileptr):
        """ Function to write self to a file such that it can be read later

        Parameters
        ----------
        fileptr : file object
            A file object to write this instance to

        """
        # Build fields to write
        # First set the sequence
        fields = [self.primer, self.diagnostic, self.reverse]

        # Now add the tags as ';' separated values
        if len(self.labels) != 0:
            fields.append(self._labelsToString())

        # Finally, write to file
        print(*fields, sep=',', file=fileptr)


class ConservedEndAmplicons:
    """
    Helper class to print sequences which have conserved ends and are all of
    equal length. Sequences are added using the add method which takes as
    arguments the forward primer, the diagnostic region, the reverse primer,
    and a label (typically a filename). Alignments can then be printed by
    calling the str() method or printing directly.
    """
    # Class level variables
    ENABLE_DOT = False
    P3_ARGS = {}

    def __init__(self, ingroup=None):
        """ Store data as a list of amplicons """
        self.amplicons = []
        self.ingroup = None
        if ingroup is not None:
            self.ingroup = frozenset(ingroup)
        self.p3 = None

    def labels(self):
        """ Returns the set of labels contained in this alignment """
        labels = set()
        for ampl in self.amplicons:
            labels |= set(ampl.labels)
        return labels

    def primerPair(self):
        """ Return the forward and reverse primers as a tuple """
        if len(self.amplicons) > 0:
            forward = self.amplicons[0].primer
            reverse = self.amplicons[0].reverse
            return (forward, reverse)
        else:
            raise ValueError("No amplicons added yet")

    def primerLength(self):
        """ Return the length of the primer region """
        if len(self.amplicons) > 0:
            return len(self.amplicons[0].primer)
        else:
            raise ValueError("No amplicons added yet")

    def diagnosticLength(self):
        """ Return the length of the diagnostic region """
        if len(self.amplicons) > 0:
            return len(self.amplicons[0].diagnostic)
        else:
            raise ValueError("No amplicons added yet")

    def ampliconLength(self):
        """ Return the length of the entire amplicon """
        if len(self.amplicons) > 0:
            return len(self.amplicons[0].sequence)
        else:
            raise ValueError("No amplicons added yet")

    def canAdd(self, other):
        """ Returns true if other can be added to this

        Parameters
        ----------
        other : Amplicon or ConservedEndAmplicons
            The instance to test addition

        Returns
        -------
        bool
            True if it can be added, else False

        Raises
        ------
        TypeError
            If input is not an Amplicon or ConservedEndAmplicons instance


        """
        # Split on type
        if isinstance(other, Amplicon):
            # Check if empty
            if len(self.amplicons) == 0:
                return True

            # Check if primers are conserved
            primer = other.primer
            reverse = other.reverse
            return (primer, reverse) == self.primerPair()
        elif isinstance(other, ConservedEndAmplicons):
            # Check if empty
            if len(self.amplicons) == 0:
                return True

            # Check if primers are conserved
            return self.primerPair() == other.primerPair()
        else:
            raise TypeError("Unexpected type input")

    def add(self, other):
        """ Add other to this instance

        Parameters
        ----------
        other : Amplicon or ConservedEndAmplicons
            The instance to add

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If input is not an Amplicon or ConservedEndAmplicons instance

        """
        # Split on input type
        if isinstance(other, Amplicon):
            # Check if this can be merged with other amplicons
            for i in range(len(self.amplicons)):
                if other == self.amplicons[i]:
                    self.amplicons[i] += other
                    break
            else:
                # Can't merge so just add to list
                self.amplicons.append(other)
        elif isinstance(other, ConservedEndAmplicons):
            # Add all amplicons to self
            for amplicon in other.amplicons:
                self.add(amplicon)
        else:
            raise TypeError("Unexpected input type")

    def diagnosticColumns(self):
        """ Returns a list of diagnostic positions """
        # Get the diagnostic sequence for each amplicon
        diagnostics = [ampl.diagnostic for ampl in self.amplicons]

        # Determine diagnostic positions
        diags = []
        for i, bases in enumerate(zip(*diagnostics)):
            if len(set(bases)) > 1:
                diags.append(i)
        return diags

    def ingroupUniqueColumns(self):
        """ Returns a list of columns unique to the ingroup """
        # Return empty list if ingroup is not set
        if self.ingroup is None:
            return []

        # Get the diagnostic sequence for each amplicon
        ingroup_diag = []
        outgroup_diag = []
        for amplicon in self.amplicons:
            # Iterate through each label
            for label in amplicon.labels:
                if label in self.ingroup:
                    ingroup_diag.append(amplicon.diagnostic)
                else:
                    outgroup_diag.append(amplicon.diagnostic)

        # Iterate through all positions
        diags = []
        for i in range(self.diagnosticLength()):
            # Get bases for ingroup and outgroup
            ingroup_bases = set([d[i] for d in ingroup_diag])
            outgroup_bases = set([d[i] for d in outgroup_diag])
            # Check if conserved in ingroup and not in outgroup
            if ingroup_bases.isdisjoint(outgroup_bases):
                diags.append(i)
        return diags

    def makeBracket(self):
        """ Make graphical bracket to output with alignment """
        # Get start and end positions of diagnostic sequences
        start = self.primerLength()
        end = start + self.diagnosticLength()

        # Create an empty bracket
        bracket = ' ' * (start-1) + '{' + '-' * (end-start) + '}'

        # Add * at diagnostic positions
        bracket = list(bracket)
        for d in self.diagnosticColumns():
            bracket[start + d] = "*"

        # Add # at unique positions in ingroup
        for d in self.ingroupUniqueColumns():
            bracket[start + d] = "#"
        return ''.join(bracket)

    def setIngroup(self, grouping):
        """ Set the 'ingroup' for the alignment, only changes str format """
        if grouping is not None:
            self.ingroup = frozenset(grouping)

    def ingroup_consensus(self):
        return self.consensus(self.ingroup)

    def consensus(self, labels=None):
        if labels is None:
            amps = self.amplicons
        else:
            amps = [x for x in self.amplicons if set(x.labels).issubset(labels)]
        forward_seq = collapse_to_iupac([x.primer for x in amps])
        diag_seq = collapse_to_iupac([x.diagnostic for x in amps])
        reverse_seq = collapse_to_iupac([x.reverse for x in amps])
        return {'forward' : forward_seq, 'diagnostic' : diag_seq, 'reverse' : reverse_seq}

    def find_primers(self):
        template = "".join(self.ingroup_consensus().values())
        self.p3 = run_primer3(template, target_start=self.primerLength(),
                              target_len=self.diagnosticLength(), **ConservedEndAmplicons.P3_ARGS)
        return self.p3['PRIMER_PAIR_NUM_RETURNED'] != 0

    def _render_primer3_stats(self):
        if self.p3 is None:
            raise ValueError('Primer3 not run before results are to be rendered.')

        left_stats = {k[14:]: v for k, v in self.p3.items() if 'PRIMER_LEFT_0_' in k}
        right_stats = {k[15:]: v for k, v in self.p3.items() if 'PRIMER_RIGHT_0_' in k}
        pair_stats = {k[14:]: v for k, v in self.p3.items() if 'PRIMER_PAIR_0_' in k}

        def render_col_names(names):
            return [x.title().replace('_', ' ') for x in names]

        def render_col_values(names):
            return [str(round(x, 5)) if isinstance(x, float) else x for x in names]

        primer_table = PrettyTable(['Direction'] + render_col_names(left_stats.keys()))
        primer_table.add_row(['Forward'] + render_col_values(left_stats.values()))
        primer_table.add_row(['Reverse'] + render_col_values(right_stats.values()))
        primer_table.align = 'l'

        pair_table = PrettyTable(render_col_names(pair_stats.keys()))
        pair_table.add_row(render_col_values(pair_stats.values()))
        pair_table.align = 'l'

        #set_trace(term_size=(80, 60))

        output = '\nPrimer statistics:\n' + \
                 primer_table.get_string(border=False) + \
                 '\n\nPair statistics:\n' + \
                 pair_table.get_string(border=False)
        return output


    def render_alignment(self):
        """ Return a string representation of an alignment """
        # Try splitting amplicons based on ingroup, outgroup
        result = []
        if self.ingroup is not None:
            in_result = []
            out_result = []
            for ampl in sorted(self.amplicons, key=lambda x: x.labels):
                if set(ampl.labels) & set(self.ingroup):
                    in_result.append(str(ampl))
                else:
                    out_result.append(str(ampl))
            result = in_result + out_result
        else:
            for ampl in sorted(self.amplicons, key=lambda x: x.labels):
                result.append(str(ampl))

        # Switch to dot alignment if enabled
        if ConservedEndAmplicons.ENABLE_DOT:
            # Mark conserved nucleotides as '.'
            top_seq = result[0]
            new_result = [top_seq]
            for seq in result[1:]:
                seq = list(seq)
                for i in range(self.ampliconLength()):
                    b0 = top_seq[i]
                    b1 = seq[i]
                    if b0 == b1:
                        seq[i] = '.'
                new_result.append(''.join(seq))
            result = new_result
        else:
            # Add diagnostic bracket to result
            result.append(self.makeBracket())

        # Add primer 3 primer annotations
        if self.p3 is not None:
            forward_seq = self.p3['PRIMER_LEFT_0_SEQUENCE']
            reverse_seq = self.p3['PRIMER_RIGHT_0_SEQUENCE']
            forward_start = self.p3['PRIMER_LEFT_0'][0]
            reverse_start = self.p3['PRIMER_RIGHT_0'][0] - self.p3['PRIMER_RIGHT_0'][1]
            forward_annot = '└' + 'Forward'.center(len(forward_seq) - 2, '─') + '┘'
            reverse_annot = '└' + 'Reverse'.center(len(reverse_seq) - 2, '─') + '┘'
            text_out = ' ' * (forward_start) + \
                       forward_annot + \
                       ' ' * (reverse_start - forward_start - len(forward_seq) + 1) + \
                       reverse_annot
            if ConservedEndAmplicons.ENABLE_DOT:
                result.append(text_out)
            else:
                result[-1] = result[-1].ljust(len(text_out))
                result[-1] = "".join([annot if bracket == ' ' else bracket
                                      for bracket, annot in zip(result[-1], text_out)])


        # Add primer3 statistics for primers
        if self.p3 is not None:
            result.append(self._render_primer3_stats())

        # Add newline to separate inputs
        result[-1] += '\n'

        # Join result with newlines and return
        return '\n'.join(result)

    def render_csv(self, sep=','):
        if len(self.amplicons) == 1:
            output = list(self.consensus().values())
        else:
            output = list(self.ingroup_consensus().values())
        if self.p3 is not None:
            output.extend(_format_p3_output(self.p3).values())
        output = [str(x) for x in output]
        return sep.join(output)

    def __str__(self):
        return self.render_alignment()

    def __len__(self):
        """ Return number of alignments """
        return len(self.amplicons)

    def __eq__(self, other):
        """ Return True if self and other have the same primer pairs """
        return self.primerPair() == other.primerPair()

    def __lt__(self, other):
        """ Return True if self comes before other in terms of primers """
        return self.primerPair() < other.primerPair()

    def __iadd__(self, other):
        """ Add other to self and return self """
        # Add other to self and return
        self.add(other)
        return self
