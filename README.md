# krisp

A python package for designing diagnostic CRISPR and other PCR-based diagnostic assays using whole genome data.

![PyPI](https://img.shields.io/pypi/v/KRISP?label=pypi%20krisp)
![PyPI - Downloads](https://img.shields.io/pypi/dm/KRISP)


## Overview

CRISPR diagnostic assays, such as SHERLOCK and DETECTR, require the design of a CRISPR guide RNA (crRNA) which contains a spacer sequence unique to the taxon of interest.
Most assays also require the design of primers for nucleic acid pre-amplification, using techniques such as PCR, RPA, or LAMP.
The increasing abundance of whole-genome sequence data can be leveraged to predict effective primers and diagnostic regions.
This data typically takes the form of whole-genome sequences (FASTA files) or variant data from reads mapped to a reference (VCF files).
The `krisp` package contains python functions and command line programs to infer diagnostic regions from either FASTA files or a VCF file plus a reference.
The `krisp_fasta` command uses a k-mer based approach on FASTA input and the `krisp_vcf` uses a sliding-window-based analysis of variants in a VCF file.
Both commands use optimized algorithms that take advantage of multiple CPUs and use minimal RAM.

`krisp_fasta` searches for diagnostic regions by converting input genome data into sorted kmer tables which are then intersected to find conserved and unique regions.
`krisp_fasta` allows specification of the desired length of the diagnostic region as well as the lengths of the conserved regions both upstream and downstream of the diagnostic regions.
This can be used to directly search for a diagnostic spacer sequence, e.g. a 28-nucleotide sequence which is conserved at every position except a single diagnostic position.
Alternatively, one can search for conserved primer regions flanking a much larger diagnostic region, e.g. 60 nucleotide diagnostic regions flanked by 32 nucleotide conserved regions upstream and downstream.
When conserved primer regions are included in the search, primers can be analyzed with Primer3 and results filtered by the presence of optimal primers. 

`krisp_vcf` scans through a VCF file with the associated reference FASTA file to find clusters of one or more diagnostic variants.
Data can be scanned to distinguish two or more groups of samples in a single execution of the command.
A sliding-window analysis of the variants searches for locations where one or more variants unique and conserved to a group of samples are surrounded by conserved sequence or variants.
Once such locations are identified, the conserved consensus sequence for the target group is inferred by applying the variants to the reference sequence.
Primer3 can then be used to identify primers in the conserved flanking regions and the regions are sorted by the predicted quality of primers.


## Requirements

`krisp` requires python version 3.6 or greater and a UNIX-like operating system which provides the `sort` command line utility.
Most operating systems have `sort` installed by default as it is required by the POSIX standard.
To see if `sort` is already installed, run:

`sort --help`

If a help menu is displayed, then `sort` is installed.
Various python packages are required as well, including `pysam`, `Bio`, `primer3-py`, and `nltk`.
These should be installed automatically when `krisp` is installed.

## Installation

### Using pip

`krisp_fasta` can be installed directly from the PyPI repository:

`pip3 install krisp`


### Downloading from github

Alternatively, `krisp_fasta` can be installed from source utilizing the Makefile in this repository.
Simply run the following command in the directory with the downloaded source files:

`make install` 

Similarly, `krisp` can be uninstalled by running:

`make uninstall`

or

`pip3 uninstall krisp`


### What is included

The `krisp` package is written in python, and when installed, provides stand-alone command line utilities called `krisp_fasta` and `krisp_vcf`.
All the functionality of `krisp_fasta` is provided through the command line and knowledge of python is not required.
In all the examples below, `krisp_fasta` is executed through the command line.

In addition to `krisp_fasta` and `krisp_vcf`, installation also provides a second command line utility called `kstream`.
`kstream` is a separate python program which extracts and parses kmers from an input FASTA file.
`krisp_fasta` utilizes `kstream` internally.
Users of `krisp_fasta` are not required to use `kstream` directly.
See examples below on `kstream` usage.

All commands have help menus which can be activated by passing the `-h` or `--help` option to the program.


## `krisp_fasta`

After the `krisp` python package has been installed, `krisp_fasta` should be available on the command line.
The help menu for `krisp_fasta` can be accessed with the `--help` option on the command line:

```
krisp_fasta --help
```

### Example 1: Searching for a spacer sequence

In example 1, we will be directly searching for a spacer sequence for CAS13a LwaCas13a, which accepts a 28-nucleotide spacer sequence and has the highest accuracy when the diagnostic SNP is at the 3rd rightmost position.
In this case, we are looking for a spacer which has the following form:

{25 conserved nts}{1 diagnostic nt}{2 conserved nts}

And we want the diagnostic nucleotide to be conserved in our "ingroup" and different in our "outgroup".

The general form of this command is as follows:

```
krisp_fasta {ingroup files} --outgroup {outgroup files} --conserved-left 25 --conserved-right 2 --diagnostic 1
```

If `krisp_fasta` was installed via the GitHub repository, then test files can be found in the "test_data/krisp_fasta" directory.
Within the "test_data/krisp_fasta" directory, we can run:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1
```

The output you see should look something like:

```
left_seq,diag_seq,right_seq
CGACAAGATACTCTCGCAGCTTGGT,M,AG
TGACGCAGATCATCCCGCGCTTACT,K,A
```

This is in the CSV format, which can be easily read by spreadsheet programs and scripting languages like R for further processing.
Instead of printing to standard output, the output can be saved to a file using the `--out_csv` option or using a bash redirect `>`:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1 --out_csv test_out.csv
```

or 

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1 > test_out.csv
```

When verbose is enabled using the `--verbose` option, the files being used and progress updates are printed during runtime:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1 --cores 4 --verbose > test_out.csv
```

In this example, all output is written to the terminal with the verbose information preceding the alignments.
Internally, the verbose data is written to stderr, whereas the alignments are written to stdout, so `test_out.csv` will not contain the progress messages.

A more human-readable alignment format can be output as well using the `--out_align` option:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1 --out_align test_align.txt
```

This produces a file with the following contents: 

```
CGACAAGATACTCTCGCAGCTTGGTCAG : ingroup0
CGACAAGATACTCTCGCAGCTTGGTAAG : ingroup1
CGACAAGATACTCTCGCAGCTTGGTGAG : outgroup0;outgroup1;outgroup2
                        {#}

TGACGCAGATCATCCCGCGCTTACTGAC : ingroup0
TGACGCAGATCATCCCGCGCTTACTTAC : ingroup1
TGACGCAGATCATCCCGCGCTTACTCAC : outgroup0;outgroup1;outgroup2
                        {#}
```

We see that 2 potential spacer sequences were found.
Each alignment is displayed as a reference sequence, followed by the alternative sequences with conserved nucleotides replaced with '.' if the `--dot-alignment` option is used:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 25 --conserved-right 2 --diagnostic 1 --dot-alignment --out_align test_align.txt
```

which produces:

```
CGACAAGATACTCTCGCAGCTTGGTCAG : ingroup0
.........................A.. : ingroup1
.........................G.. : outgroup0;outgroup1;outgroup2

TGACGCAGATCATCCCGCGCTTACTGAC : ingroup0
.........................T.. : ingroup1
.........................C.. : outgroup0;outgroup1;outgroup2
```

To the right of each sequence is the filename(s) in which the sequence was found, with multiple filenames implying that all of these files contain the same sequence.
Note that it is possible to have multiple sequences listed for a single file, which may be common for repetitive regions or diploid (or higher) organisms.


### Example 2: Searching for conserved primer regions

`krisp_fasta` can search for conserved primer regions which span both the ingroup and outgroup and can use Primer3 to design primers.
The general formula is the same as in example 1, except that our conserved and diagnostic region will be longer.

Let us assume that we are looking for an amplicon which is 100 nucleotides in total length, and contains at least 30 conserved nucleotides on each end where we can design primers.
To find these regions, we can run the following `krisp_fasta` command:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved-left 30 --conserved-right 30 --diagnostic 40 --dot-alignment --primer3 --out_align test_align.txt
```

This command was designed to be deliberately similar to that in example 1. 
For convenience, `krisp_fasta`, when the conserved flanking regions are the same, we can simply specify `--conserved 30`.
We can also specify the total amplicon length as opposed to calculating the diagnostic length: `--amplicon 100`.
In which case, the above command is equivalent to:

```
krisp_fasta ingroup*.gz --outgroup outgroup*.gz --conserved 30 --amplicon 100 --dot-alignment --primer3 --out_align test_align.txt
```

Both commands should produce output like the following in "test_align.txt":

```
ACGCACAAGGACAAGTGCCACTAAACCAGCCAGCCCTGACGCAGATCATCCCGCGCTTACTGACCAAGCTGCGAGAGTATCTTGTCGATGGGAACGATAG : ingroup0
.............................................................T...................................... : ingroup1
.............................................................C...................................... : outgroup0;outgroup1;outgroup2
   └────────Forward─────────┘                                           └────────Reverse────────┘

Primer statistics:
 Direction  Penalty  Sequence                    Tm        Gc Percent  Self Any Th  Self End Th  Hairpin Th  Position Penalty  End Stability  Template Mispriming       Template Mispriming Th   
 Forward    7.74706  CACAAGGACAAGTGCCACTAAACCAG  64.24706  50.0        0.0          2.14676      0.0         0.0               4.0            -1.7976931348623157e+308  -1.7976931348623157e+308 
 Reverse    6.43757  TCGTTCCCATCGACAAGATACTCTC   61.93757  48.0        0.0          0.0          37.5163     0.0               3.2            -1.7976931348623157e+308  -1.7976931348623157e+308 

Pair statistics:
 Penalty   Compl Any Th  Compl End Th  Product Size  Product Tm  Product Tm Oligo Tm Diff  T Opt A   Template Mispriming      
 14.18463  0.0           0.0           94            84.32116    22.38359                  62.70608  -1.7976931348623157e+308 
```

In this example, `krisp_fasta` finds two regions which specify our input constraints, although only the first is shown in the output above.
Note that these alignments are actually reverse complements. 
`krisp_fasta` assumes that the input genomes are double-stranded and creates the corresponding kmers for both strands, thus finding two solutions.

The addition of the `--primer3` option causes Primer3 to be run for each potential region.
When `--primer3` is used, only regions for which it was possible to find primers are included in the output.
Additionally, the Primer3 statistics are included in the alignment and CSV output.

### Example 3: Using `krisp_fasta` in creative ways

Users are encouraged to use `krisp_fasta` in any way that they seem fit.
In this example, we will be using krisp to find all 60-mers which are conserved across a group of input files.
To do this, we will specify a conserved region of 30-nts in each direction and a diagnostic region of 0-nts.

```                                                                             
krisp_fasta ingroup*.gz outgroup*.gz --conserved 30 --diagnostic 0                                        
```

The result of which is:

```
left_seq,diag_seq,right_seq
ACGCACAAGGACAAGTGCCACTAAACCAGC,,CAGCCCTGACGCAGATCATCCCGCGCTTAC
AGTAAGCGCGGGATGATCTGCGTCAGGGCT,,GGCTGGTTTAGTGGCACTTGTCCTTGTGCG
CGCACAAGGACAAGTGCCACTAAACCAGCC,,AGCCCTGACGCAGATCATCCCGCGCTTACT
GTAAGCGCGGGATGATCTGCGTCAGGGCTG,,GCTGGTTTAGTGGCACTTGTCCTTGTGCGT
```

In this case, there is actually no distinction between the ingroup and outgroup files as these are only required to differ in the diagnostic region which we specified to be of 0 length.


### Example 4: Using `kstream` to extract kmers

`kstream` is a flexible and powerful tool for extracting and parsing kmers from an input fasta file or input stream of sequences.
Some available options are listed below, but see the help menu for a complete list.

`--kmers` allows the user to specify the length (or lengths) or kmers to extract.

`--disallow` allows the user to omit any kmers containing disallowed nucleotides. See `--allow` for its inverse.

`--map-softmask` will map all masked nucleotides (lowercase characters) back to uppercase characters in all kmers.

`--omit-softmask` will omit any kmer containing soft-masked nucleotides.

`--complements` will add the complementary kmers to the stream.

`--canonicals` will only print out the "canonical" kmers, e.g. those which come first alphabetically.

`--sort` will sort the resulting stream.

`--parallel` will utilize multiple cores to speed up extraction.


#### Extracting complete sequences

By default, `kstream` will parse a fasta file and print out every sequence if no kmer length is provided. For example:

`kstream test_DNA.fasta`

Will print every sequence entry in the input file and disregard all header names (lines starting with '>').
If an input file is not given, then `kstream` listens to the standard input pipe.
The previous command can also be written as:

`cat test_DNA.fasta | kstream`


#### Sorted 6-mers, without 'N' or 'n', removing softmasking, in canonical form

This example utilizes a few of the common command line arguments: For example,
to extract all 6-mers which don't contain "N" or "n", in sorted order:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort`

If we only wanted the "canonical" kmers, then we could add the `--canonicals` flag:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort --canonicals`

or pass the previous stream into a new `kstream` instance:

`kstream test_DNA.fasta --kmers 6 --disallow "Nn" --sort | kstream --canonicals`

Both are equivalent but the former is likely to be faster when cpus are
limiting.


## `krisp_vcf` 

After the `krisp` python package has been installed, `krisp_vcf` should be available on the command line.
The help menu for `krisp_vcf` can be accessed with the `--help` option on the command line:

```
krisp_vcf --help
```

### Basic usage

`krisp_vcf` requires VCF input (compressed or not) and the reference genome in FASTA format used to create it.
It also requires a metadata CSV file with at least two columns: one containing sample IDs that match the column names in the VCF file and one with arbitrary values corresponding to which group each sample belongs to.
The first few rows of the metadata file for the example data looks like this (extra spaces added to make it easier to read):

```
sample_id,group
9D1,NA1
BS2014-584,NA1
BS96,NA1
CC1008,EU1
CC1011,EU1
CC1033,EU1
```

Files used in the examples below can be found in the "test_data/krisp_vcf" of this repository.
You can use them to test `krisp_vcf` using the example commands below.
Here is a minimal command:

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1`

This will scan the variants in "variants.vcf" and look for diagnostic regions that distinguish samples labeled as "NA1", "NA2", or "EU1" from all other samples.
The VCF file is provided with the option `--vcf` since it could also be piped from standard input.
Which group each sample is assigned to is defined by the "metadata.csv" file.
Since no output file was specified, CSV output is printed to the standard output and status updates are printed to standard error.
When run on the command line, both the output (standard out) and the status updates (standard error) are printed to the screen, which is not ideal most of the time.
We can save the CSV output to a file by redirecting the output using `>` or using the `--out_csv` option:

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv`

or

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 > test_out.csv`

While the program is running it prints status updates like this:

```
Undiagnostic| Unconserved | No primers  | NA1         | NA2         | EU1         
20050       | 91          | 14          | 6           | 22          | 17          
```

These numbers are counts of variant clusters considered and why they were determined to not be diagnostic or which group they are diagnostic for.

The output of the program is saved in "test_out.csv", which has many columns:
* "group": The group the region is diagnostic for
* "chrom": The chromosome on the reference database 
* "n_diag": The number of bases that are diagnostic in the region
* A group of columns that contain the reference positions of primers, the diagnostic regions, and nearby sequences
* A group of columns containing inferred consensus sequence for the target group for primers, the diagnostic regions, and nearby sequences
* A group of columns containing the output of Primer3

This format is ideal for post-processing with bioinformatic tools and programming languages like R.

### The alignment output

A more human-readable output format can also be generated by supplying the `--out_align` option: 

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv --out_align test_align.txt`

This prints each output as a plain text alignment with the Primer3 statistics:

```
## Phyram_PR-102_s0001:179343-179465 is diagnostic for EU1

Reference : ctgggatactactaccatttgaccaattcaaggccatgtatgcctacgaccgctactagcgttggaaagttg  T  cgcttgacttgtgccacagcaaatg  T  tgacatttgga  G  aattggtctccg
      EU1 : ........................................................................<C13>.........................<C12>...........<A12>............
      NA1 : ........................................................................  T6 .........................  T6 ...........  G6 ............
      NA2 : ........................................................................ T10 ......................... T10 ...........  G8 ............
   oligos : ctgggatactactaccatttgaccaattc                  gaccgctactagcgttggaaagttg  C  cg                         C  tgacatttgga  A  aattggtctccg
            └─────── Left primer ───────┘                  └─────────── crRNA ────────────┘                       └───────── Right primer ────────┘

Primer statistics:
 Direction  Penalty  Sequence                       Tm        Gc Percent  Self Any Th  Self End Th  Hairpin Th  Position Penalty  End Stability  Template Mispriming       Template Mispriming Th   
 Forward    2.04051  ctgggatactactaccatttgaccaattc  61.54051  41.37931    0.0          0.0          31.61709    0.0               2.17           -1.7976931348623157e+308  -1.7976931348623157e+308 
 Reverse    5.34751  cggagaccaattTtccaaatgtcaG      60.84751  44.0        0.0          0.0          0.0         0.0               3.51           -1.7976931348623157e+308  -1.7976931348623157e+308 

Pair statistics:
 Penalty  Compl Any Th  Compl End Th  Product Size  Product Tm  Product Tm Oligo Tm Diff  T Opt A   Template Mispriming      
 7.38802  5.70965       4.50833       123           83.24808    22.40056                  61.62791  -1.7976931348623157e+308 
```

This alignment output contains the reference sequence and the inferred consensus sequence for each group being distinguished.
The alignment uses the following syntax:

* Dots (`.`): positions in the group consensus sequence that match the reference
* Lowercase letters: Sequence derived from the reference FASTA rather than variants from the VCF
* Uppercase letters in reference/primers: variants from the VCF, potentially conserved in all groups
* Series of uppercase letters with numbers (e.g. `T6`): A variant with the number of samples supporting that variant. For example, `T6AG2` would mean 6 samples had a `T` and 2 samples had a `AG` insertion at this position.
* Series of uppercase letters with numbers flanked by angle brackets (e.g. `<C12>`): Variants diagnostic for the group they are printed in and the number of samples that have this variant. 


### Primer3 options

Most of the commonly used settings of Primer3 can be set on the command line:

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv --out_align test_align.txt --amp_size 50 100 --gc_clamp 3`

Note how the addition of `--amp_size 50 100 --gc_clamp 3` changes the size of amplicons returned and the kind of primers selected:

```
## Phyram_PR-102_s0001:182872-182966 is diagnostic for NA2

Reference : caatatacctcctcaaacatcacaccctcAacttcatttattgtcggCccaacgcatgG  T  CcgcttcaatccgcatactaagtatttcCaaccg A  
      EU1 : ........................................................... T14 ..................................G14 
      NA1 : ...........................................................  T6 ..................................<A6>
      NA2 : ...........................................................<G12>..................................G12 
   oligos : caatatacctcctcaaacatcacaccc       catttattgtcggcccaacgcatgg  G  cc        ccgcatactaagtatttccaaccg g  
            └────── Left primer ──────┘       └─────────── crRNA ────────────┘        └────── Right primer ──────┘
```

Other Primer3 settings that can be changed include primer length, annealing temperature, GC content, maximum temperature of secondary structures (e.g. primer dimers).
See the help menu with `krisp_vcf --help` for more Primer3 options. 

### Input subsetting and parallel processing

`krisp_vcf` can be made to run much faster by subsetting the input to specific chromosomes or positions of interest.
To subset to a specific reference position (for all chromosomes), use the `--pos` option: 

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv --pos 100000 200000`

To subset to a specific chromosome (a single sequence in the reference FASTA file, identified by its header), use the `--chroms` option: 

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv --chroms Phyram_PR-102_s0001`

The `--pos` and `--chroms` options can be used together to subset to a specific position on a specific chromosome.

If multiple processors are available on the computer running `krisp_vcf`, parallel processing can be enabled by specifying the number of processors to use using the `--cores` option:

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 NA2 EU1 --out_csv test_out.csv --cores 4`


### Quality filtering

VCF files can contain erroneous variants depending on how they have been quality filtered.
`krisp_vcf` can do its own quality filtering of variants, allowing unfiltered variants to be used as input.
By default, variants must occur in 5 samples per group, be represented by at least 10 reads in each of those samples, and have a genotype quality score of 40.
These defaults can be seen in the help menu displayed by `krisp_vcf --help`.
Custom values can be set using the `--min_samples`, `--min_reads`, and `--min_geno_qual` options: 

`krisp_vcf metadata.csv reference.fasta --vcf variants.vcf --groups NA1 EU1 --out_csv test_out.csv --min_samples 3 --min_reads 30`

### License

This work is subject to the MIT License.

### Feedback and contributions

We would like to hear about users’ thoughts on the package and any errors they run into. Please report errors, questions or suggestions on the issues tab in this github repository. We also welcome contributions via a Github pull request. 
