#!/usr/bin/env python

from pprint import pprint

import optparse
from Bio import SeqIO
import numpy as np
from la import larry

def detect_fastq_quality_encoding(filename, nnuc = 50000):
    '''Given a file name, parse the file as a FASTQ file and attempt
    to determine the quality encoding used. Returns either "sanger",
    "solexa", or "illumina".

    By default, only the first 50,000 nucleotides will be considered.'''
    # I am aware that this method is very un-general. I don't know how
    # to make it more general.

    # The following are used to eliminate various possibilities.
    # https://secure.wikimedia.org/wikipedia/en/wiki/Fastq#Encoding
    # offset_ranges = { sanger: (33,73), solexa: (59,104), illumina: (64,104) }
    sanger_offset = 33
    sanger_max_threshold = 40
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
        if 'sanger' in possible_encodings and max_seen > sanger_max_threshold:
            possible_encodings.remove('sanger')
        if 'solexa' in possible_encodings and min_seen < solexa_min_threshold:
            possible_encodings.remove('solexa')
        if 'illumina' in possible_encodings and min_seen < illumina_min_threshold:
            possible_encodings.remove('illumina')
        if len(possible_encodings) == 1:
            return possible_encodings.pop()
        elif len(possible_encodings) == 0:
            raise Exception("Could not identify FASTQ file: eliminated all possible encodings.")
        nuc_count += len(record)
        if nuc_count > nnuc:
            break
    # If no Illumina-encoded quality less than zero has been seen,
    # then eliminate solexa and return illumina.
    if possible_encodings == set('solexa', 'illumina'):
        return 'illumina'
    else:
        return 'solexa'

def fastq_parse_autodetect_format(filename):
    '''Same as SeqIO.parse(filename, "fastq"), except that instead of
    "fastq", the encoding format will be automatically detected.'''
    encoding_format = "fastq-%s" % (detect_fastq_quality_encoding(filename),)
    return SeqIO.parse(filename, encoding_format)

def seq_to_data_triplets(seq_record):
    '''Takes a sequence and returns a list of tuples (BASE, POSITION, QUALITY).

    Position is ZERO-based, so you can use it as an array index.'''
    return zip(seq_record.seq,
               xrange(0:len(seq_record)),
               seq_record.letter_annotations["phred_quality"])

def quality_historgram_by_position_and_base(seqio):
    '''Returns a big 3-D labeled array.

    Assumes that all sequences are the same length, and that the only
    ambiguous base is N.'''
    bases = ( 'A', 'T', 'C', 'G', 'N', )
    max_quality = 40
    # Get sequence length from first sequence.
    first_seq = seqio.next()
    read_length = len(first_seq)

    hist = np.zeros((len(bases),read_length, max_quality))


    first_triplets = seq_to_data_triplets(first_seq)

    read_length = 100





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

if __name__ == "__main__":
    (options, args) = parse_options()
    pprint((options, args,))
