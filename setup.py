import os
from setuptools import setup, find_packages

def read_file(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "fastqident",
    version = "0.1dev",
    packages = find_packages(),
    entry_points = {
        'console_scripts': [
            'fastqident = fastqident.main:plac_call_main',
        ],
    },

    # metadata for upload to PyPI
    author = "Ryan C. Thompson",
    author_email = "rct@thompsonclan.org",
    description = "Guess the quality encoding system used by FASTQ sequence files.",
    license = "BSD",
    keywords = ("biopython", "bioinformatics", "fastq", "sequencing"),
    #url = "http://example.com/HelloWorld/",   # project home page, if any
    long_description=read_file('README.md'),
    requires = ('plac', 'placsupport'),
)
