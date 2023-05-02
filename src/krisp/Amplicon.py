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
    def __init__(self, ingroup=None):
        """ Store data as a list of amplicons """
        self.amplicons = []
        self.ingroup = None
        if ingroup is not None:
            self.ingroup = frozenset(ingroup)

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

    def __str__(self):
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

        # Join result with newlines and return
        return '\n'.join(result)

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
