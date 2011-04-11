from __future__ import print_function

try:
    import plac
    from placsupport import print_argument_error
    from placsupport.types import *
except ImportError, e:
    print("You must install the plac and placsupport modules to use the fastqident command-line script.")
    raise e

# Entry point
def plac_call_main():
    return plac.call(main)

@plac.annotations(
    # (helptext, kind, abbrev, type, choices, metavar)
    fastq=('The FASTQ files to identify', 'option'),
    max_quality=('Assumed maximum possible quality value', 'option', 'm', positive_int),
    nnuc=('Number of nuelcotides to sample from each file', 'option', 'n', positive_int),
    initskip=('How many sequences to skip at the beginning of each file', 'option', 'i', nonneg_int),
    skip=('How many sequences to skip between each sampled sequence', 'option', 's', nonneg_int),
    possible_encodings=('Comma-separated list of possible quality encodings', 'option', 'e'),
    sanger_min=('Minimum ASCII value of Sanger-encoded qualities', 'option', 'g', positive_int),
    solexa_min=('Minimum ASCII value of Solexa-encoded qualities', 'option', 'x', positive_int),
    illumina_min=('Minimum ASCII value of Illumina-encoded qualities', 'option', 'l', positive_int),
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
    """Detect the quality encoding of FASTQ sequence files."""
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
