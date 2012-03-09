# fastqident: Guess the quality encoding system used by FASTQ sequence files.

The FASTQ sequence format stores biological sequences, along with
sequencing quality values for each element in each sequence.
Unfortunately, different sources encode the quality values in
different ways. The different encodings are sparsely documented, and
telling them apart is
[confusing](http://en.wikipedia.org/wiki/FASTQ_format#Encoding),
especially because many files are techincally valid in multiple
encodings. However, in practice it is generally possible to make a
good guess as to which encoding was used based on the observed range
of ASCII character. The purpose of this module is to do this guessing
for you. It should get the answer right as long as the question isn't
too hard.

## Installation

fastqident is distributed as a standard Python package. You can download the tarball here: https://github.com/DarwinAwardWinner/fastqident/tarball/master

Simply install it your normal Python package installer. Probably `pip`
or `easy_install` or setup.py.

## Usage

### From the command line

This package includes a script by the same name that takes a list of
fastq files on the command line and tries to identify the encoding of
each. Usage:

    $ fastqident read1.fastq read2.fastq
    {'read1.fastq': 'illumina', 'read2.fastq': 'illumina'}

### From python code

Import the constructor:

    from fastqident import FastqQualityIdentifier

Create an identifier with default values:

    id_default = FastqQualityIdentifier()

Create an identifier with custom options:

    id_custom = FastqQualityIdentifier(max_quality=50, nnuc=5000, start=50, skip=50)

The identifier class supplies three methods: `detect_encoding`,
`detect_encoding_safe`, and `detect_encodings`. Here is some example
code using them:

    # Identify a single fastq file (returns a string)
    filename = "read1.fastq"
    file_encoding = id_default.detect_encoding(filename)
    print "%s has quality encoding %s" % (filename, file_encodings)

    # Identify a list of files (returns a dict)
    filenames = ("read1.fastq", "read2.fastq", "read3.fastq")
    file_encodings = id_custom.detect_encodings(filenames)
    for fname, fenc in file_encodings:
        print "%s has quality encoding %s" % (fname, fenc)

## Errata

fastqident assumes that the sequences near the start of the file will
provide a reasonable sample of the range of quality values in the
whole file. If this assumption fails, then fastqident may fail to
correctly identify things. You can increase `-n`, `-i`, and `-s` from
their default values to take a larger sample across a greater fraction
of the fastq file, at the cost of taking longer.

Examples of misidentified files are welcome.
