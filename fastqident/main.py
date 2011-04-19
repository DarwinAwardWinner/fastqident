# This contains the main function for the command-line script, along
# with its associated argument parsing infrastructure.

from __future__ import print_function

try:
    import plac
    from placsupport import argument_error
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
    max_quality=('Upper bound on the quality values you expect to see in the file. The tighter you can make this bound, the more reliable the identification will be. However, if you set this bound too low, you will get misidentifications.', 'option', 'm', positive_int),
    nnuc=('Number of nucleotides to sample from each file', 'option', 'n', positive_int),
    initskip=('How many sequences to skip at the beginning of each file. You may wish to use this because the first sequences in an Illumina data set are often uninformative.', 'option', 'i', nonneg_int),
    skip=('How many sequences to skip between each sampled sequence. This allows you to sample from a larger range of the fastq file.', 'option', 's', nonneg_int),
    possible_encodings=('Comma-separated list of possible quality encodings. By default, all three possibilities are checked, but if you have already ruled out one of them, simply list the other two.', 'option', 'e'),
    sanger_min=('Minimum ASCII value of Sanger-encoded qualities', 'option', 'g', positive_int),
    solexa_min=('Minimum ASCII value of Solexa-encoded qualities', 'option', 'x', positive_int),
    illumina_min=('Minimum ASCII value of Illumina-encoded qualities', 'option', 'l', positive_int),
    allow_empty_file_list=('Do not produce an error if no fastq files are passed on the command line. This is potentially useful for running in a pipeline.', 'flag','z'),
    )
def main(nnuc = 50000,
         initskip = 0,
         skip = 4,
         possible_encodings = 'sanger,solexa,illumina',
         sanger_min = 33,
         solexa_min = 59,
         illumina_min = 64,
         max_quality = 40,
         allow_empty_file_list=False,
         *fastq):
    """Detect the quality encoding of FASTQ sequence files."""
    if not allow_empty_file_list and len(fastq) < 1:
        argument_error('Need at least one fastq file to operate on')
    possible_encodings = set(map(str.strip, possible_encodings.split(",")))
    known_encodings = set(('sanger', 'solexa', 'illumina'))
    unknown_encodings = possible_encodings.difference(known_encodings)
    if unknown_encodings:
        argument_error('The only known encodings are sanger, solexa, and illumina. You supplied the following unknown encodings:\n\n%s' % (",".join(sorted(unknown_encodings))))
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
