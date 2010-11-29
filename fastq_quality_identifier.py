#!/usr/bin/env python

# TODO: convert this file into a fqident package
# TODO: Rewrite to be more efficient using "from Bio.SeqIO.QualityIO import FastqGeneralIterator"

from __future__ import print_function
from Bio import SeqIO
from copy import copy

class FastqQualityIdentifier(object):
    def __init__(self, nnuc = 50000, max_quality = 40,
                 possible_encodings = set(('sanger', 'solexa', 'illumina')),
                 sanger_min = 33, solexa_min = 59, illumina_min = 64):
        # For info on the default values, see:
        # https://secure.wikimedia.org/wikipedia/en/wiki/Fastq#Encoding
        self.nnuc = nnuc
        if not self.nnuc or self.nnuc <= 0:
            self.nnuc = False
        self.max_quality = max_quality
        self.possible_encodings = set(possible_encodings)
        self.sanger_min = sanger_min
        self.solexa_min = solexa_min
        self.illumina_min = illumina_min
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
        # Any valid illumina-format fastq is also valid sanger, so parsing
        # as sanger format will always succeed.
        for record in SeqIO.parse(filename, "fastq-sanger"):
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
    import sys
    # Use pprint if available. Otherwise, use regular print
    try:
        from pprint import pprint
        printfunc = pprint
    except ImportError:
        printfunc = print
    x = FastqQualityIdentifier()
    printfunc(x.detect_encodings(sys.argv[1:]))
