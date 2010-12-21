#!/usr/bin/env python

# TODO: convert this file into a fqident package
# TODO: Rewrite to be more efficient using "from Bio.SeqIO.QualityIO import FastqGeneralIterator"

from __future__ import print_function
from Bio import SeqIO
from copy import copy
from itertools import islice

class FastqQualityIdentifier(object):
    def __init__(self, max_quality = 40, nnuc = 50000, start = 0, skip = 4,
                 possible_encodings = set(('sanger', 'solexa', 'illumina')),
                 sanger_min = 33, solexa_min = 59, illumina_min = 64):
        # For info on the default values, see:
        # https://secure.wikimedia.org/wikipedia/en/wiki/Fastq#Encoding
        self.nnuc = nnuc
        if not self.nnuc or self.nnuc <= 0:
            self.nnuc = None
        self.max_quality = max_quality
        self.possible_encodings = set(possible_encodings)
        self.sanger_min = sanger_min
        self.solexa_min = solexa_min
        self.illumina_min = illumina_min
        self.start = start
        self.skip = skip
        # Computed values
        self.solexa_threshold = self.solexa_min - self.sanger_min
        self.illumina_threshold = self.illumina_min - self.sanger_min

    def detect_encoding(self, filename):
        '''Given a file name, parse the file as a FASTQ file and
        attempt to determine the quality encoding. Returns either
        "sanger", "solexa", or "illumina".

        max_quality is the assumed maximum quality value that any
        nucleotide can have. The default is 40, the anecdotal maximum
        quality for an Illumina sequencing dataset.

        nnuc is the number of nucleotides to read. It is assumed that
        in reading this many nucleotides starting at the beginning of
        the file, the script will encounter the full range of possible
        quality values, from which it will then make a descision. The
        default is 50_000.

        This docstring is for a function that was subsequently
        rewritten into a class. The documentation is obsolete.'''
        # I am aware that this method is very un-general. I don't know how
        # to make it more general.
        possible_encodings = copy(self.possible_encodings)
        # Min starts high and works it way down. Vice versa for max.
        min_seen = 128
        max_seen = 0
        nuc_count = 0
        seq_count = 0
        # Any valid illumina-format fastq is also valid sanger, so parsing
        # as sanger format will always succeed.
        seqio = SeqIO.parse(filename, "fastq-sanger")
        seqio_slice = islice(seqio, self.start, None, self.skip)
        for record in seqio_slice:
            seq_count += 1
            qualities = record.letter_annotations["phred_quality"]
            min_seen = min(min_seen, min(qualities))
            max_seen = max(max_seen, max(qualities))
            # Eliminate possibilities
            if 'sanger' in possible_encodings and max_seen > self.max_quality:
                possible_encodings.remove('sanger')
            if 'solexa' in possible_encodings and min_seen < self.solexa_threshold:
                possible_encodings.remove('solexa')
            if 'illumina' in possible_encodings and min_seen < self.illumina_threshold:
                possible_encodings.remove('illumina')
            # Check if we finished early
            if len(possible_encodings) == 1:
                return possible_encodings.pop()
            elif len(possible_encodings) == 0:
                raise Exception("Could not identify FASTQ file %s: eliminated all possible encodings." % (filename,))
            if self.nnuc:
                nuc_count += len(record)
                if nuc_count > self.nnuc:
                    break
        # If no Illumina-encoded quality less than zero has been seen,
        # then eliminate solexa and return illumina.
        if min_seen >= self.illumina_threshold:
            return 'illumina'
        else:
            return 'solexa'

    def _detect_encoding_safe(self, filename):
        '''Same as detect_encoding, but returns None on IOError.'''
        try:
            return self.detect_encoding(filename)
        except IOError:
            return None

    def detect_encodings(self, filenames):
        '''Detect the quality encodings of each of a list of files.

        Returns a dict with filenames as keys and encoding styles as
        values.'''
        return dict([ (f, self._detect_encoding_safe(f)) for f in filenames ])

if __name__ == "__main__":
    # import sys
    import plac
    def print_help(fun):
        """Print help text inferred by plac for fun"""
        plac.parser_from(fun).print_help()

    def print_argument_error(fun, message):
        """Print error message, followed by usage message for fun"""
        print("\nError: %s\n" % (message, ))
        print_help(fun)
    @plac.annotations(
        # (helptext, kind, abbrev, type, choices, metavar)
        fastq=('The FASTQ files to identify', 'option'),
        max_quality=('Assumed maximum possible quality value', 'option', 'm'),
        nnuc=('Number of nuelcotides to sample from each file', 'option', 'n'),
        initskip=('How many sequences to skip at the beginning of each file', 'option', 'i'),
        skip=('How many sequences to skip between each sampled sequence', 'option', 's'),
        possible_encodings=('Comma-separated list of possible quality encodings', 'option', 'e'),
        sanger_min=('Minimum ASCII value of Sanger-encoded qualities', 'option', 'g'),
        solexa_min=('Minimum ASCII value of Solexa-encoded qualities', 'option', 'x'),
        illumina_min=('Minimum ASCII value of Illumina-encoded qualities', 'option', 'l'),
        allow_empty_file_list=('Do not prooduce an error for zero fastq files. This is potentially useful for running in a pipeline.', 'flag','z'),
        )
    def main(max_quality = 40,
             nnuc = 50000,
             initskip = 0,
             skip = 4,
             possible_encodings = 'sanger,solexa,illumina',
             sanger_min = 33,
             solexa_min = 59,
             illumina_min = 64,
             allow_empty_file_list=False,
             *fastq):
        if not allow_empty_file_list and len(fastq) < 1:
            print_argument_error(main, 'Need at least one fastq file to operate on')
        possible_encodings = set(map(str.strip, possible_encodings.split(",")))
        known_encodings = set(('sanger', 'solexa', 'illumina'))
        unknown_encodings = possible_encodings.difference(known_encodings)
        if unknown_encodings:
            print_argument_error(main, 'The only known encodings are sanger, solexa, and illumina. You supplied the following unknown encodings:\n\n%s' % (",".join(sorted(unknown_encodings))))
        if type(sanger_min) == str:
            sanger_min = ord(sanger_min)
        if type(solexa_min) == str:
            solexa_min = ord(solexa_min)
        if type(illumina_min) == str:
            illumina_min = ord(illumina_min)
        x = FastqQualityIdentifier(max_quality, nnuc, initskip, skip,
                                   possible_encodings, sanger_min,
                                   solexa_min, illumina_min)
        # Use pprint if available. Otherwise, use regular print
        try:
            from pprint import pprint
            printfunc = pprint
        except ImportError:
            printfunc = print
        printfunc(x.detect_encodings(fastq))
    plac.call(main)
