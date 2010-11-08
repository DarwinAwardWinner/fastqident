#!/usr/bin/env python

from pprint import pprint

import optparse
from Bio import SeqIO

def parse_options():
    parser = optparse.OptionParser(usage='%prog -o OUTFILE FASTQ_FILES...', version='%prog 0.1')
    parser.add_option('-o', '--output-file',
                      help='The name of the output file to write. This will be a numpy data file. The default is to take the first input file, strip ".fastq" from the end, and add ".npz".',
                      action='store',
                      dest='outfile',
                      default=False,)
    parser.add_option('-q', '--quality-encoding',
                      help='The method used for quality encoding in the input files. This can be "sanger", "solexa", or "illumina". You can also specify "auto", in which case the script will do its best to auto-detect the encoding used. This will take a bit of extra time per file. The default is "auto".',
                      action='store',
                      dest='qual_enc',,
                      default='auto')
    (options, args) = parser.parse_args()
    return (options, args)

def detect_quality_encoding(filename):
    pass

if __name__ == "__main__":
    (options, args) = parse_options()
    pprint((options, args,))
