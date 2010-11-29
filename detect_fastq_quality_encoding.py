#!/usr/bin/env python

from pprint import pprint

import sys
from Bio import SeqIO

def detect_fastq_quality_encoding(filename, max_quality = 40, nnuc = 50000):
    '''Given a file name, parse the file as a FASTQ file and attempt
    to determine the quality encoding. Returns either "sanger",
    "solexa", or "illumina".

    max_quality is the assumed maximum quality value that any
    nucleotide can have. The default is 40, the anecdotal maximum
    quality for an Illumina sequencing dataset.

    nnuc is the number of nucleotides to read. It is assumed that in
    reading this many nucleotides starting at the beginning of the
    file, the script will encounter the full range of possible quality
    values, from which it will then make a descision. The default is
    50_000.'''
    # I am aware that this method is very un-general. I don't know how
    # to make it more general.

    # The following are used to eliminate various possibilities.
    # https://secure.wikimedia.org/wikipedia/en/wiki/Fastq#Encoding
    # offset_ranges = { sanger: (33,73), solexa: (59,104), illumina: (64,104) }
    sanger_offset = 33
    max_quality = 40
    solexa_min_threshold = 59 - sanger_offset
    illumina_min_threshold = 64 - sanger_offset

    possible_encodings = set(('sanger', 'solexa', 'illumina'))

    # Min starts high and works it way down. Vice versa for max
    min_seen = 128 - sanger_offset
    max_seen = 0
    nuc_count = 0
    # Any valid illumina-format fastq is also valid sanger, so parsing
    # as sanger format will always work.
    for record in SeqIO.parse(filename, "fastq-sanger"):
        qualities = record.letter_annotations["phred_quality"]
        min_seen = min(min_seen, min(qualities))
        max_seen = max(max_seen, max(qualities))
        if 'sanger' in possible_encodings and max_seen > max_quality:
            possible_encodings.remove('sanger')
        if 'solexa' in possible_encodings and min_seen < solexa_min_threshold:
            possible_encodings.remove('solexa')
        if 'illumina' in possible_encodings and min_seen < illumina_min_threshold:
            possible_encodings.remove('illumina')

        if len(possible_encodings) == 1:
            print "Min: %s; Max: %s" % (min_seen + 33, max_seen + 33)
            return possible_encodings.pop()
        elif len(possible_encodings) == 0:
            raise Exception("Could not identify FASTQ file: eliminated all possible encodings.")
        if nnuc > 0:
            nuc_count += len(record)
            if nuc_count > nnuc:
                break

    print "Min: %s; Max: %s" % (min_seen + 33, max_seen + 33)
    # If no Illumina-encoded quality less than zero has been seen,
    # then eliminate solexa and return illumina.
    if min_seen >= illumina_min_threshold:
        return 'illumina'
    else:
        return 'solexa'

if __name__ == "__main__":
    print "fastq-%s" % (detect_fastq_quality_encoding(sys.argv[1], nnuc=5000),)
