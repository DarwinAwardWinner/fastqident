#!/usr/bin/env python

from pprint import pprint

# import optparse
from Bio import SeqIO
# import numpy as np
# from la import larry
from fastq_quality_identifier import FastqQualityIdentifier
from itertools import chain,izip

def fastq_parse_autodetect(filename):
    '''Same as SeqIO.parse(filename, "fastq"), except that instead of
    "fastq", the encoding format will be automatically detected.'''
    encoding_format = "fastq-%s" % (FastqQualityIdentifier().detect_encoding(filename),)
    return SeqIO.parse(filename, encoding_format)

def seqrecord_to_data_triplets(seq_record):
    '''Takes a sequence and returns a list of tuples (BASE, POSITION, QUALITY).

    Position is ZERO-based, so you can use it as an array index.'''
    return izip(seq_record.seq,
                xrange(0,len(seq_record)),
                seq_record.letter_annotations["phred_quality"])

def seqio_to_data_triplets(seqio, limit = 1000):
    '''Takes a SeqIO object and returns tuples of (BASE, POS, QUAL).'''
    count = 0
    for rec in seqio:
        print "Reading", rec.id
        seq = rec.seq
        qual = rec.letter_annotations["phred_quality"]
        for i in xrange(0,len(seq)):
            count += 1
            yield (seq[i], i, qual[i])
            if limit and count >= limit:
                raise StopIteration()

    # return chain.from_iterable([ seqrecord_to_data_triplets(seq) for seq in seqio ])

def quality_historgram_by_position_and_base(seqio, limit = 1000):
    hist = {}
    triplet_iter = seqio_to_data_triplets(seqio, limit=limit)
    if limit:
        count = 0
    for t in triplet_iter:
        try:
            hist[t] += 1
        except KeyError:
            hist[t] = 1
        if limit:
            count += 1
            if count >= limit:
                break
    return hist




# def parse_options():
#     parser = optparse.OptionParser(usage='%prog -o OUTFILE FASTQ_FILES...', version='%prog 0.1')
#     parser.add_option('-o', '--output-file',
#                       help='The name of the output file to write. This will be a numpy data file. The default is to take the first input file, strip ".fastq" from the end, and add ".npz".',
#                       action='store',
#                       dest='outfile',
#                       default=False,)
#     parser.add_option('-q', '--quality-encoding',
#                       help='The method used for quality encoding in the input files. This can be "sanger", "solexa", or "illumina". You can also specify "auto", in which case the script will do its best to auto-detect the encoding used. This will take a bit of extra time per file. The default is "auto".',
#                       action='store',
#                       dest='qual_enc',
#                       default='auto')
#     (options, args) = parser.parse_args()
#     return (options, args)

if __name__ == "__main__":
    import sys
    seqio = fastq_parse_autodetect(sys.argv[1])
    for x in seqio_to_data_triplets(seqio, limit=10000):
        print x
    # pprint(quality_historgram_by_position_and_base(seqio), limit=500)
