# krisp
A lightweight kmer-based algorithm for designing diagnostic CRISPR assays
using genome data.

## Overview
CRISPR diagnostic assays, such as SHERLOCK and DETECTR, require the designing
of a crRNA which contains a spacer sequence unique to the pathogen in question.
*krisp* searches for these unique spacers by converting input genome data into
sorted kmer tables which are then intersected to find conserved and unique
regions.

*krisp* allows the user to specify the desired length of the diagnostic
region as well as the lengths of the conserved regions both upstream and
downstream of the diagnostic region. This can be used to directly search for a
diagnostic spacer sequence, e.g. a 28-nucleotide sequence which is conserved
at every position except a single diagnostic position. Alternatives, one can
search for conserved primer regions flanking a much larger diagnostic region,
e.g. a 60 nucleotide diagnostic regions flanked by 32 nucleotide conserved
regions upstream and downstream.

## Requirements
*krisp* requires python version 3.6 or greater and a UNIX-like operating
system which provides the `sort` command line utility. Most operating
systems have `sort` installed by default as it is required by the POSIX
standard. To see if `sort` is already installed, run:

`sort --help`

## Installation
### Using pip
*krisp* can be installed directly from the python repository
(testpypi for now):

`pip3 install -i https://test.pypi.org/simple/krisp`

### Downloading from github
Alternatively, *krisp* can be installed from source utilizing the
provided Makefile. Simply run in the root directory:

`make install` 

Similarly, krisp can be uninstalled by running:

`make uninstall`

or

`pip3 uninstall krisp`

### What is included
The *krisp* package is written in python, and when installed, provides a
stand-alone command line utility also called *krisp*. All of the functionality
of *krisp* is provided through the command line and knowledge of
python is not required. In all the examples below, *krisp* is executed through
the command line.

In addition to *krisp*, installation also provides a second command line 
utility called *kstream*. *kstream* is a seperate python program which extracts
and parses kmers from an input fasta file. *krisp* utilizes *kstream*
internally and is provided as a courtesy. Users of *krisp* are never
required to utilize *kstream* in any *krisp* workflow. See example 4 below the
*kstream* usage.

Both *krisp* and *kstream* provide help menus which can be activated by passing
the `-h` or `--help` option to the program.

## Examples
### Example 1: Searching for a spacer sequence
In example 1, we will be directly searching for a spacer sequence for
LwaCas13a which accepts a 28-nucleotide spacer sequence and has the highest
accuracy when the diagnostic SNP is at the 3rd rightmost position.
In this case, we are looking for a spacer which has the following form:

{25 conserved nts}{1 diagnostic nt}{2 conserved nts}

And we want the diagnostic nucleotide to be conserved in our "ingroup" and 
different in our "outgroup".

The general form of this command is as follows:

```
krisp {ingroup files} \
        --outgroup {outgroup files} \
        --conserved-left 25 \
        --conserved-right 2 \
        --diagnostic 1 \
        --verbose \
        --parallel {number of processors}

```

If *krisp* was installed via the github repository, then test files can be
found in the testdata/ directory. Within the testdata/ directory, we can run:

```
krisp ingroup*.gz \
        --outgroup outgroup*.gz \
        --conserved-left 25 \
        --conserved-right 2 \
        --diagnostic 1 \
        --verbose \
        --parallel 4

```

The output you see should look something like:

```
verbose
Finding kmer-based diagnostic regions for:
(0) ingroup0.fasta.gz
(1) ingroup1.fasta.gz
With this as an outgroup:
(0) outgroup0.fasta.gz
(1) outgroup1.fasta.gz
(2) outgroup2.fasta.gz

Extracting 28-mers from ingroup0.fasta.gz and saving to /tmp/tmpx2q99a9h/ingroup0.28mers
Extracting 28-mers from ingroup1.fasta.gz and saving to /tmp/tmpx2q99a9h/ingroup1.28mers
Extracting 28-mers from outgroup0.fasta.gz and saving to /tmp/tmpx2q99a9h/outgroup0.28mers
Extracting 28-mers from outgroup1.fasta.gz and saving to /tmp/tmpx2q99a9h/outgroup1.28mers
Extracting 28-mers from outgroup2.fasta.gz and saving to /tmp/tmpx2q99a9h/outgroup2.28mers
=> Extracted and sorted 19,660 28-kmers from ingroup1.fasta.gz in 0.07 seconds
=> Extracted and sorted 19,660 28-kmers from outgroup1.fasta.gz in 0.07 seconds
=> Extracted and sorted 19,660 28-kmers from outgroup0.fasta.gz in 0.07 seconds
=> Extracted and sorted 19,660 28-kmers from ingroup0.fasta.gz in 0.11 seconds
=> Extracted and sorted 19,660 28-kmers from outgroup2.fasta.gz in 0.11 seconds
Merging: /tmp/tmpx2q99a9h/outgroup2.28mers + /tmp/tmpx2q99a9h/outgroup1.28mers -> /tmp/tmpx2q99a9h/tmpmzv_vxqx using 1 cores
Merging: /tmp/tmpx2q99a9h/outgroup0.28mers + /tmp/tmpx2q99a9h/ingroup1.28mers -> /tmp/tmpx2q99a9h/tmpw8bkaq0h using 1 cores
=> Merged /tmp/tmpx2q99a9h/outgroup2.28mers + /tmp/tmpx2q99a9h/outgroup1.28mers -> /tmp/tmpx2q99a9h/tmpmzv_vxqx in 0.30 seconds
=> Merged /tmp/tmpx2q99a9h/outgroup0.28mers + /tmp/tmpx2q99a9h/ingroup1.28mers -> /tmp/tmpx2q99a9h/tmpw8bkaq0h in 0.31 seconds
Merging: /tmp/tmpx2q99a9h/ingroup0.28mers + /tmp/tmpx2q99a9h/tmpw8bkaq0h -> /tmp/tmpx2q99a9h/tmpamgmksp9 using 1 cores
=> Merged /tmp/tmpx2q99a9h/ingroup0.28mers + /tmp/tmpx2q99a9h/tmpw8bkaq0h -> /tmp/tmpx2q99a9h/tmpamgmksp9 in 0.15 seconds
Merging: /tmp/tmpx2q99a9h/tmpmzv_vxqx + /tmp/tmpx2q99a9h/tmpamgmksp9 -> /tmp/tmpx2q99a9h/tmpp9dicqvd using 1 cores
=> Merged /tmp/tmpx2q99a9h/tmpmzv_vxqx + /tmp/tmpx2q99a9h/tmpamgmksp9 -> /tmp/tmpx2q99a9h/tmpp9dicqvd in 0.01 seconds

Building alignments ... 
> Alignment 0
CGACAAGATACTCTCGCAGCTTGGTCAG : ingroup0
.........................A.. : ingroup1
.........................G.. : outgroup0;outgroup1;outgroup2
> Alignment 1
TGACGCAGATCATCCCGCGCTTACTGAC : ingroup0
.........................T.. : ingroup1
.........................C.. : outgroup0;outgroup1;outgroup2
=> Found 2 alignments in 0.60 seconds

```

When verbose is enabled, *krisp* will list the files being used as the ingroup
and outgroup, and print progress information during runtime. Note that
the verbose output will likely differ from that above as temporary files with
random names are generated to store intermediate results. These tempory files
are stored in the computers designated tempory directory, often `/tmp`, and are
automatically cleaned up after the program exits.

In the alignment section, we see that 2 potential spacer sequences were found,
labeled as Alignment 0 and Alignment 1 respectively. Each alignment is
displayed as a reference sequence, followed by the alternative sequences with
conserved nucleotides replaced with '.'. To the right of each sequence is the
filename(s) in which the sequence was found. With multiple filenames implying
that all of these files contain the same sequence. Note that it is possible to
have multiple sequences listed for a single file, which may be common for
repetitive regions or diploid (or higher) organisms.

### Example 2: Searching for conserved primer regions
In example 2, we will be using *krisp* to search for conserved primer regions
which span both the ingroup and outgroup. The general formula is the same as
in example 1, except that our conserved and diagnostic region will be
longer.

Let us assume that we are looking for an amplicon which is 100 nucleotides in
total length, and contains at least 30 conserved nucleotides on each end where
we will design primers. To find these regions, we can run the following *krisp*
command:

```
krisp ingroup*.gz \
        --outgroup outgroup*.gz \
        --conserved-left 30 \
        --conserved-right 30 \
        --diagnostic 40 \
        --verbose \
        --parallel 4
```

This command was designed to be deliberately similar to that in example 1. For
convenience, *krisp* has many other command line options. For instance, when
the conserved flanking regions are the same, one can simply specify
`--conserved 30` or `-c 30`. We can also specify the total amplicon length as
opposed to calculating the diagnostic length: `--amplicon 100` or `-a 100`. In
which case, the above command is equivalent to:

```
krisp ingroup*.gz \
        --outgroup outgroup*.gz \
        --conserved 30 \
        --amplicon 100 \
        --verbose \
        --parallel 4
```

The output of which should look something like:

```
Finding kmer-based diagnostic regions for:
(0) ingroup0.fasta.gz
(1) ingroup1.fasta.gz
With this as an outgroup:
(0) outgroup0.fasta.gz
(1) outgroup1.fasta.gz
(2) outgroup2.fasta.gz

Extracting 100-mers from ingroup0.fasta.gz and saving to /tmp/tmpxabqwtjd/ingroup0.100mers
Extracting 100-mers from ingroup1.fasta.gz and saving to /tmp/tmpxabqwtjd/ingroup1.100mers
Extracting 100-mers from outgroup0.fasta.gz and saving to /tmp/tmpxabqwtjd/outgroup0.100mers
Extracting 100-mers from outgroup1.fasta.gz and saving to /tmp/tmpxabqwtjd/outgroup1.100mers
Extracting 100-mers from outgroup2.fasta.gz and saving to /tmp/tmpxabqwtjd/outgroup2.100mers
=> Extracted and sorted 18,220 100-kmers from outgroup1.fasta.gz in 0.13 seconds
=> Extracted and sorted 18,220 100-kmers from outgroup2.fasta.gz in 0.13 seconds
=> Extracted and sorted 18,220 100-kmers from ingroup0.fasta.gz in 0.13 seconds
=> Extracted and sorted 18,220 100-kmers from outgroup0.fasta.gz in 0.13 seconds
=> Extracted and sorted 18,220 100-kmers from ingroup1.fasta.gz in 0.13 seconds
Merging: /tmp/tmpxabqwtjd/outgroup2.100mers + /tmp/tmpxabqwtjd/outgroup1.100mers -> /tmp/tmpxabqwtjd/tmppdmq3pvq using 1 cores
Merging: /tmp/tmpxabqwtjd/outgroup0.100mers + /tmp/tmpxabqwtjd/ingroup1.100mers -> /tmp/tmpxabqwtjd/tmpwft9lb3x using 1 cores
=> Merged /tmp/tmpxabqwtjd/outgroup2.100mers + /tmp/tmpxabqwtjd/outgroup1.100mers -> /tmp/tmpxabqwtjd/tmppdmq3pvq in 0.28 seconds
=> Merged /tmp/tmpxabqwtjd/outgroup0.100mers + /tmp/tmpxabqwtjd/ingroup1.100mers -> /tmp/tmpxabqwtjd/tmpwft9lb3x in 0.28 seconds
Merging: /tmp/tmpxabqwtjd/ingroup0.100mers + /tmp/tmpxabqwtjd/tmpwft9lb3x -> /tmp/tmpxabqwtjd/tmpu2cg3j3t using 1 cores
=> Merged /tmp/tmpxabqwtjd/ingroup0.100mers + /tmp/tmpxabqwtjd/tmpwft9lb3x -> /tmp/tmpxabqwtjd/tmpu2cg3j3t in 0.14 seconds
Merging: /tmp/tmpxabqwtjd/tmppdmq3pvq + /tmp/tmpxabqwtjd/tmpu2cg3j3t -> /tmp/tmpxabqwtjd/tmpqh8stw2c using 1 cores
=> Merged /tmp/tmpxabqwtjd/tmppdmq3pvq + /tmp/tmpxabqwtjd/tmpu2cg3j3t -> /tmp/tmpxabqwtjd/tmpqh8stw2c in 0.02 seconds

Building alignments ... 
> Alignment 0
ACGCACAAGGACAAGTGCCACTAAACCAGCCAGCCCTGACGCAGATCATCCCGCGCTTACTGACCAAGCTGCGAGAGTATCTTGTCGATGGGAACGATAG : ingroup0
.............................................................T...................................... : ingroup1
.............................................................C...................................... : outgroup0;outgroup1;outgroup2
> Alignment 1
CTATCGTTCCCATCGACAAGATACTCTCGCAGCTTGGTCAGTAAGCGCGGGATGATCTGCGTCAGGGCTGGCTGGTTTAGTGGCACTTGTCCTTGTGCGT : ingroup0
......................................A............................................................. : ingroup1
......................................G............................................................. : outgroup0;outgroup1;outgroup2
=> Found 2 alignments in 0.59 seconds

```

In this example, *krisp* finds two regions which specify our input constraints.
Note that these alignments are actually reverse complements. *krisp* assumes
that the input genomes are double-stranded and creates the corresponding kmers
for both strands, thus finding two solutions.

### Example 3: Using *krisp* in creative ways
Users are encouraged to use *krisp* in any way that they seem fit. In this
example, we will be using krisp to find all 60-mers which are critically 
conserved across a group of input files. To do this, we will specify a
conserved region of 30-nts in each direction and a diagnostic region of 0-nts.

```                                                                             
krisp ingroup*.gz \                                                             
        --outgroup outgroup*.gz \                                               
        --conserved 30 \                                                        
        --diagnostic 0 \                                                        
        --verbose \                                                             
        --parallel 4                                                           
```

The result of which is:

```
Finding kmer-based diagnostic regions for:
(0) ingroup0.fasta.gz
(1) ingroup1.fasta.gz
With this as an outgroup:
(0) outgroup0.fasta.gz
(1) outgroup1.fasta.gz
(2) outgroup2.fasta.gz

Extracting 60-mers from ingroup0.fasta.gz and saving to /tmp/tmpbx3fbnt8/ingroup0.60mers
Extracting 60-mers from ingroup1.fasta.gz and saving to /tmp/tmpbx3fbnt8/ingroup1.60mers
Extracting 60-mers from outgroup0.fasta.gz and saving to /tmp/tmpbx3fbnt8/outgroup0.60mers
Extracting 60-mers from outgroup1.fasta.gz and saving to /tmp/tmpbx3fbnt8/outgroup1.60mers
Extracting 60-mers from outgroup2.fasta.gz and saving to /tmp/tmpbx3fbnt8/outgroup2.60mers
=> Extracted and sorted 19,020 60-kmers from outgroup0.fasta.gz in 0.09 seconds
=> Extracted and sorted 19,020 60-kmers from outgroup2.fasta.gz in 0.09 seconds
=> Extracted and sorted 19,020 60-kmers from outgroup1.fasta.gz in 0.09 seconds
=> Extracted and sorted 19,020 60-kmers from ingroup0.fasta.gz in 0.13 seconds
=> Extracted and sorted 19,020 60-kmers from ingroup1.fasta.gz in 0.13 seconds
Merging: /tmp/tmpbx3fbnt8/outgroup2.60mers + /tmp/tmpbx3fbnt8/outgroup1.60mers -> /tmp/tmpbx3fbnt8/tmpvi16uh4y using 1 cores
Merging: /tmp/tmpbx3fbnt8/outgroup0.60mers + /tmp/tmpbx3fbnt8/ingroup1.60mers -> /tmp/tmpbx3fbnt8/tmpmrtyeubg using 1 cores
=> Merged /tmp/tmpbx3fbnt8/outgroup2.60mers + /tmp/tmpbx3fbnt8/outgroup1.60mers -> /tmp/tmpbx3fbnt8/tmpvi16uh4y in 0.31 seconds
=> Merged /tmp/tmpbx3fbnt8/outgroup0.60mers + /tmp/tmpbx3fbnt8/ingroup1.60mers -> /tmp/tmpbx3fbnt8/tmpmrtyeubg in 0.31 seconds
Merging: /tmp/tmpbx3fbnt8/ingroup0.60mers + /tmp/tmpbx3fbnt8/tmpmrtyeubg -> /tmp/tmpbx3fbnt8/tmprw0a31pb using 1 cores
=> Merged /tmp/tmpbx3fbnt8/ingroup0.60mers + /tmp/tmpbx3fbnt8/tmpmrtyeubg -> /tmp/tmpbx3fbnt8/tmprw0a31pb in 0.14 seconds
Merging: /tmp/tmpbx3fbnt8/tmpvi16uh4y + /tmp/tmpbx3fbnt8/tmprw0a31pb -> /tmp/tmpbx3fbnt8/tmpqk8ur264 using 1 cores
=> Merged /tmp/tmpbx3fbnt8/tmpvi16uh4y + /tmp/tmpbx3fbnt8/tmprw0a31pb -> /tmp/tmpbx3fbnt8/tmpqk8ur264 in 0.01 seconds

Building alignments ... 
> Alignment 0
ACGCACAAGGACAAGTGCCACTAAACCAGCCAGCCCTGACGCAGATCATCCCGCGCTTAC : ingroup0;ingroup1;outgroup0;outgroup1;outgroup2
> Alignment 1
AGTAAGCGCGGGATGATCTGCGTCAGGGCTGGCTGGTTTAGTGGCACTTGTCCTTGTGCG : ingroup0;ingroup1;outgroup0;outgroup1;outgroup2
> Alignment 2
CGCACAAGGACAAGTGCCACTAAACCAGCCAGCCCTGACGCAGATCATCCCGCGCTTACT : ingroup0;ingroup1;outgroup0;outgroup1;outgroup2
> Alignment 3
GTAAGCGCGGGATGATCTGCGTCAGGGCTGGCTGGTTTAGTGGCACTTGTCCTTGTGCGT : ingroup0;ingroup1;outgroup0;outgroup1;outgroup2
=> Found 4 alignments in 0.61 seconds

```

In this case, there is actually no distinction between the ingroup and outgroup
files as these are only required to differ in the diagnostic region which we 
specified to be of 0 length.

### Example 4: Using *kstream* to extract kmers
*kstream* is an extreamly flexible and powerful tool for extracting and parsing
kmers from an input fasta file or input stream of sequences. Some of the
available options are listed below, but see the help menu for a complete list.


`--kmers` allows the user to specify the length (or lengths)
or kmers to extract.

`--disallow` allows the user to omit any kmers containing
disallowed nucleotides. See `--allow` for its inverse.

`--map-softmask` will map all masked nucleotides (lowercase characters) back to
uppercase characters in all kmers.

`--omit-softmask` will omit any kmer containing soft-masked nucleotides

`--complements` will add the complementary kmers to the stream

`--canonicals` will only print out the "canonical" kmers, e.g. those which come
first alphabetically.

`--sort` will sort the resulting stream

`--parallel` will utilize multiple cores to speed up extraction 

#### Extracting complete sequences
By default, *kstream* will parse a fasta file and print out every sequence if
no kmer length is provided. For example:

`kstream test_DNA.fasta`

Will print out every sequence entry in the input file and disregard all header
names (lines starting with '>'). If an input file is not given, then *kstream*
listens to the standard input pipe.
The previous command can also be written as:

`cat test_DNA.fasta | kstream`

#### Sorted 6-mers, without 'N' or 'n', removing softmasking, in canonical form
This example utilizes a few of the common command line arguments: For example,
to extract all 6-mers which don't contain "N" or "n", in sorted order:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort`

If we only wanted the "canonical" kmers, then we could add the `--canonicals`
flag:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort --canonicals`

or pass the previous stream into a new *kstream* instance:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort | kstream --canonicals`

Both are equivalent but the former is likely to be faster when cpus are
limiting.
