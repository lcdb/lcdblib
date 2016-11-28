#!/usr/bin/env python
""" Converts between chromosome names. """
import os
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
import pkg_resources
from tempfile import NamedTemporaryFile
import gzip

import pandas as pd
import pysam
import pybedtools
from Bio import SeqIO

from lcdblib.logger import logger


def arguments():
    """ Function to pull command line arguments """

    DESCRIPTION = """\
    Converts between chromosome names.

    In D. melanogaster there are four common varieties of chromosomes names that
    are used.

    * FlyBase style: 2L, 2R, 3L, etc.
    * UCSC style: chr2L, chr2R, chr3L, etc.
    * RefSeq style: NT_033779.5, NT_033778.4, NT_037436.4, etc.
    * GenBank style: AE014134.6, AE013599.5, AE014296.5, etc.

    Each style has its benefits: FlyBase is the standard, UCSC is compatible
    with the genome browser, and RefSeq and GenBank are the most explicit. It
    is common to come across files that use these different styles of
    chromosome name, this tool aims to easily convert chromosome names in a
    variety of file format.

    """

    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     formatter_class=Raw)

    parser.add_argument("--from", dest="orig", action='store', required=True,
                        choices=['FlyBase', 'UCSC', 'RefSeq', 'GenBank'],
                        help="Current type of chromosome name.")

    parser.add_argument("--to", dest="new", action='store', required=True,
                        choices=['FlyBase', 'UCSC', 'RefSeq', 'GenBank'],
                        help="The type of chromosome name wanted.")

    parser.add_argument("--fileType", dest="type", action='store',
                        required=True,
                        choices=['SAM', 'BAM', 'BED', 'GFF', 'GTF', 'FASTA'],
                        help="What is the input format.")

    parser.add_argument("-i", "--input", dest="input", action='store',
                        required=True, help="Input file to convert.")

    parser.add_argument("-o", "--output", dest="output", action='store',
                        required=True, help="Output file.")

    parser.add_argument("--debug", dest="debug", action='store_true',
                        required=False, help="Enable debug output.")

    return parser.parse_args()


def import_conversion(f, t):
    """
    Import NCBI conversion table.

    Parameters
    ----------
    f: str {FlyBase, UCSC, NCBI}
        The current chromosome format.
    t: str {FlyBase, UCSC, NCBI}
        The desired chromosome format.

    Returns
    -------
    dict: Mapping {f: t}

    """

    # Get location of the file
    gz = pkg_resources.resource_filename('lcdblib', 'data/GCF_000001215.4.assembly.txt.gz')

    # Import file
    df = pd.read_csv(gz, compression='gzip', comment='#', sep='\t', header=None)

    # Make mapping from keyword to column number
    # col#   Header
    # 0      Sequence-Name
    # 1      Sequence-Role
    # 2      Assigned-Molecule
    # 3      Assigned-Molecule-Location/Type
    # 4      GenBank-Accn
    # 5      Relationship
    # 6      RefSeq-Accn
    # 7      Assembly-Unit
    # 8      Sequence-Length
    # 9      UCSC-style-name
    mapping = {'FlyBase': 0, 'UCSC': 9, 'GenBank': 4, 'RefSeq': 6}

    # This table has a small error were FlyBase mitochondrion_genome is MT
    df.replace('MT', 'mitochondrion_genome', inplace=True)

    return {k: v for k, v in df[[mapping[f], mapping[t]]].values}


def decompress(input):
    tmp = NamedTemporaryFile()
    with gzip.open(input, 'rb') as fh:
        tmp.file.writelines(fh.readlines())

    return tmp.name

def pysam_convert(input, output, kind, mapper):
    """
    Use pysam to convert chromosomes in BAM and SAM files.

    pysam uses a header to define chromosomes, then each read is just mapped
    back to this header. Only the header needs to be modified, and reads
    need to be written to the new output file which uses this header.

    """
    # Determine SAM or BAM flags
    if kind == 'BAM':
        flag_in = 'rb'
        flag_out = 'wb'
    elif kind == 'SAM':
        flag_in = 'r'
        flag_out = 'w'

    curr = pysam.AlignmentFile(input, flag_in)

    # Change chromosome in the header
    header = curr.header
    for chrom in header['SQ']:
        chrom['SN'] = mapper[chrom['SN']]

    with pysam.AlignmentFile(output, flag_out, header=header) as OUT:
        for read in curr:
            OUT.write(read)


def convertFeature(f, mapper):
    f.chrom = mapper[f.chrom]
    return f


def pybedtools_convert(input, output, mapper):
    """ Use pybedtools to convert chromosomes in BED, GTF, or GFF. """
    pybedtools.BedTool(input).each(convertFeature, mapper).saveas(output)


def fasta_convert(input, output, mapper):
    """ Uses Biopython.SeqIO to convert FASTA headers. """
    with open(output, 'w') as OUT:
        for seq in SeqIO.parse(input, 'fasta'):
            seq.description = seq.description.replace(seq.id, mapper[seq.id])
            seq.name = mapper[seq.name]
            seq.id = mapper[seq.id]
            SeqIO.write(seq, OUT, 'fasta')


def main():
    # Import commandline arguments.
    args = arguments()

    # Get mapping dict
    mapper = import_conversion(args.orig, args.new)

    # Extract gziped files
    if args.input.endswith('.gz'):
        input = decompress(args.input)
    else:
        input = args.input

    if (args.type == 'BAM') | (args.type == 'SAM'):
        pysam_convert(input, args.output, args.type, mapper)
    elif (args.type == 'BED') | (args.type == 'GFF') | (args.type == 'GTF'):
        pybedtools_convert(input, args.output, mapper)
    elif (args.type == 'FASTA'):
        fasta_convert(input, args.output, mapper)

