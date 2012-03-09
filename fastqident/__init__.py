#!/usr/bin/env python

# TODO: Rewrite to be more efficient using "from Bio.SeqIO.QualityIO
# import FastqGeneralIterator"

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
        '''Given a fastq file name, guess the quality encoding

        Returns either "sanger", "solexa", or "illumina".'''
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
        seqio_slice = islice(seqio, self.start, None, self.skip + 1)
        for record in seqio_slice:
            seq_count += 1
            qualities = record.letter_annotations["phred_quality"]
            min_seen = min(min_seen, min(qualities))
            max_seen = max(max_seen, max(qualities))
            # Eliminate possibilities
            if 'sanger' in possible_encodings and max_seen > self.max_quality:
                possible_encodings.remove('sanger')
            if 'solexa' in possible_encodings and min_seen < self.solexa_threshold:
                return 'sanger'
            if 'illumina' in possible_encodings and min_seen < self.illumina_threshold:
                possible_encodings.remove('illumina')
            # Check if we finished early
            if len(possible_encodings) == 1:
                return possible_encodings.pop()
            elif len(possible_encodings) == 0:
                raise ValueError("Could not identify FASTQ file %s: eliminated all possible encodings." % (filename,))
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

    def detect_encoding_safe(self, filename):
        '''Same as detect_encoding, but does not raise Exceptions.

        Instead, it can return "INVALID".'''
        try:
            return self.detect_encoding(filename)
        except ValueError, e:
            return "INVALID"

    def detect_encodings(self, filenames):
        '''Detect the quality encodings of each of a list of files.

        Returns a dict with filenames as keys and encoding styles as
        values.'''
        return dict((f, self.detect_encoding_safe(f)) for f in filenames )

__default_identifier = FastqQualityIdentifier()

detect_encoding = __default_identifier.detect_encoding
detect_encoding_safe = __default_identifier.detect_encoding_safe
detect_encodings = __default_identifier.detect_encodings
