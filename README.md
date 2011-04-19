# fastqident: Guess the quality encoding system used by FASTQ sequence files.

The FASTQ sequence format stores biological sequences, along with
sequencing quality values for each element in each sequence.
Unfortunately, different sources encode the quality values in
different ways. The purpose of this module is to make an educated
guess at which encoding style was used for a particular fastq file, so
you don't have to.

## Usage

### From the command line

This package includes a script by the same name that takes a list of
fastq files on the command line and tries to identify the encoding of
each. Usage:

    $ fastqident read1.fastq read2.fastq
    {'read1.fastq': 'illumina', 'read2.fastq': 'illumina'}

### From python code

Create an identifier with default values:

    id_default = FastqQualityIdentifier()

Create an identifier with custom options:

    id_custom = FastqQualityIdentifier(max_quality=50, nnuc=5000, start=50, skip=50)

Use the identifier:

    # Identify a single fastq file (returns a string)
    filename = "read1.fastq"
    file_encoding = id_default.detect_encoding(filename)
    print "%s has quality encoding %s" % (filename, file_encodings)

    # Identify a list of files (returns a dict)
    filenames = ("read1.fastq", "read2.fastq", "read3.fastq")
    file_encodings = id_custom.detect_encodings(filenames)
    for fname, fenc in file_encodings:
        print "%s has quality encoding %s" % (fname, fenc)

