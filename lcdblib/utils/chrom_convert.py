#!/usr/bin/env python
""" Converts between chromosome names. """
import os
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
import pkg_resources

from lcdblib.logger import logger


def arguments():
    """ Function to pull command line arguments """

    DESCRIPTION = """\
    Converts between chromosome names.

    In D. melanogaster there are three common varieties of chromosomes names that
    are used.

    * FlyBase style: 2L, 2R, 3L, etc.
    * UCSC style: chr2L, chr2R, chr3L, etc
    * NCBI style:

    Each style has its benefits: FlyBase is the standard, UCSC is compatible with
    the genome browser, NCBI is the most explicit. It is common to come across
    files that use these different styles of chromosome name, this tool aims to
    easily convert chromosome names in a variety of file format.

    """

    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     formatter_class=Raw)

    parser.add_argument("--from", dest="from", action='store', required=True,
                        choices=['FlyBase', 'UCSC', 'NCBI'],
                        help="Current type of chromosome name.")

    parser.add_argument("--to", dest="to", action='store', required=True,
                        choices=['FlyBase', 'UCSC', 'NCBI'],
                        help="The type of chromosome name wanted.")

    parser.add_argument("--fileType", dest="type", action='store',
                        required=True,
                        choices=['SAM', 'BAM', 'BED', 'BigWig'],
                        help="What is the input format.")

    parser.add_argument("-i", "--input", dest="input", action='store',
                        required=True, help="Input file to convert.")

    parser.add_argument("-o", "--output", dest="output", action='store',
                        required=True, help="Output file.")

    parser.add_argument("--debug", dest="debug", action='store_true',
                        required=False, help="Enable debug output.")

    return parser.parse_args()


def main():
    args = arguments()
    print(pkg_resources.resource_string('data', 'GCF_000001215.4.assembly.txt.gz'))
