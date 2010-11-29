#!/usr/bin/env python

from pprint import pprint

from Bio import SeqIO
# import numpy as np
# from la import larry
from fastq_quality_identifier import FastqQualityIdentifier
from itertools import chain,izip,imap

def fastq_parse_autodetect(filename):
    '''Same as SeqIO.parse(filename, "fastq"), except that instead of
    "fastq", the encoding format will be automatically detected.'''
    encoding_format = "fastq-%s" % (FastqQualityIdentifier().detect_encoding(filename),)
    return SeqIO.parse(filename, encoding_format)

def seqrecord_to_data_triplets(seq_record):
    '''Takes a sequence and returns a list of tuples (BASE, POSITION, QUALITY).

    Position is ZERO-based, so you can use it as an array index.'''
    return izip(seq_record.seq,
                xrange(1,len(seq_record)+1),
                seq_record.letter_annotations["phred_quality"])

def seqio_to_data_triplets(seqio):
    '''Takes a SeqIO object and returns tuples of (BASE, POS, QUAL).'''
    return chain.from_iterable(imap(seqrecord_to_data_triplets, seqio))

def quality_historgram_by_position_and_base(seqio, limit = 100000000):
    hist = {}
    triplet_iter = seqio_to_data_triplets(seqio)
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

def hist_to_array(hist):
    def histline_to_list(k,v):
        a = list(k)
        a.append(v)
        return a
    return [ histline_to_list(k,v) for k,v in sorted(hist.iteritems()) ]


def array_to_tsv(array):
    return "\n".join([ "\t".join([ str(s) for s in x]) for x in array ])

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

if __name__ == '__main__':
    import sys
    import plac

    @plac.annotations(
        # (helptext, kind, abbrev, type, choices, metavar)
        filename=('The input FASTQ file', 'positional', None, str, None, 'FASTQ-FILE'),
        limit=('Number of nucleotides to read (0 for unlimited)', 'option', 'l', float, None, 'COUNT'),
        outfile=('Output file (default stdout)', 'option', 'o', str, None, 'OUTPUT-FILE'),
        qenc=('Quality Encoding of the input file (default autodetect)', 'option', 'q', str, ['sanger', 'solexa', 'illumina',], 'QENC'))
    def main(filename, limit=1e9, outfile=sys.stdout, qenc=None):
        if limit <= 0:
            limit = None
        if outfile == '-':
            outfile = sys.stdout
        if outfile != sys.stdout:
            outfile = file(outfile, 'w')
        if qenc:
            seqio = SeqIO.parse(filename, 'fastq-%s' % (qenc, ))
        else:
            seqio = fastq_parse_autodetect(filename)
        hist = quality_historgram_by_position_and_base(seqio, limit)
        hist_array = [ ("base", "position", "quality", "count"), ]
        hist_array.extend(hist_to_array(hist))
        hist_tsv = array_to_tsv(hist_array)
        outfile.write(hist_tsv)

    plac.call(main)
