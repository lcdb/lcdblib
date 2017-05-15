#!/usr/bin/env python
""" Quick and dirty Fastqc parser """
import os
import pandas as pd
from io import StringIO
from zipfile import ZipFile
from tempfile import TemporaryDirectory

from lcdblib.logger import logger

logger.setLevel(10)


class FastQC(object):
    def __init__(self, id, filename):
        """ Parse a FastQC data file """
        self.filename = filename
        self.id = id
        self.blocks = {}

        self._parser()

    def _parser(self):
        with open(self.filename, 'r') as fh:
            currBlock = None
            for row in fh:
                if row:
                    if row.startswith('##FastQC'):
                        self.version = row.split()[1]
                    elif row == '>>END_MODULE\n':
                        if len(currBlock) > 1:
                            cb = FastQCBlock(currBlock)
                            self.blocks[cb.id] = cb
                        currBlock = None
                    elif row.startswith('>>'):
                        currBlock = [row.lstrip('>>').rstrip()]
                    else:
                        currBlock.append(row.rstrip())

    def __getitem__(self, key):
        return self.blocks[key]

    def keys(self):
        return self.blocks.keys()

    def values(self):
        return self.blocks.values()

    def items(self):
        return self.blocks.items()

    @classmethod
    def parse_from_zip(cls, id, filename):
        """Parse from ZipFile."""
        with TemporaryDirectory() as tmp:
            with ZipFile(filename, 'r') as archive:
                fname = [x.filename for x in archive.filelist if 'fastqc_data.txt' in x.filename][0]
                archive.extract(fname, path=tmp)
                fq = cls(id, os.path.join(tmp, fname))
        return fq


class FastQCBlock(object):
    def __init__(self, block):
        self.id = block[0].split('\t')[0]
        self.status = block[0].split('\t')[1]

        if self.id == 'Sequence Duplication Levels':
            block[2] = block[2].lstrip('#')
            block_str = '\n'.join(block[2:])
            self.df = pd.read_csv(StringIO(block_str), sep='\t', index_col=0)
        else:
            block[1] = block[1].lstrip('#')
            block_str = '\n'.join(block[1:])
            self.df = pd.read_csv(StringIO(block_str), sep='\t', index_col=0)


def split_ranges(df):
    """Split ranges into bases.

    Fastqc sometimes collapses bases into ranges, this splits them back out.
    """
    rows = []
    for i, row in df.iterrows():
        try:
            if '-' in i:
                start, end = [int(x) for x in i.split('-')]
                for j in range(start, end + 1):
                    curr_row = row.copy()
                    curr_row.name = j
                    rows.append(curr_row)
            else:
                row.name = int(i)
                rows.append(row)
        except TypeError:
            rows.append(row)

    df = pd.concat(rows, axis=1).T
    df.index.name = 'base'
    return df


def parse_fastqc(sample, file, field=''):
    """Parse fastqc zip file.

    Takes a zip file and makes a FastQC object.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    field: str
        Name of specific Fastqc section to return. Look at
        lcdblib.parse.fastqc.FastQC.keys() for a list of possible names.

    Returns
    -------
    lcdblib.parse.fastqc.FastQC or lcdblib.parse.fastqc.FastQCBlock if `field`
    is provided.

    """
    if field:
        return FastQC.parse_from_zip(sample, file)[field]
    else:
        return FastQC.parse_from_zip(sample, file)


def parse_fastqc_per_seq_quality(sample, file):
    """Parse fastqc Per Seq Quality.
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Per sequence quality scores')
    df = fqc.df
    wide = df.T
    wide['sample'] = sample
    return wide.set_index('sample')


def parse_fastqc_per_base_seq_quality(sample, file):
    """Parse fastqc Per Base Quality
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Per base sequence quality')
    df = fqc.df['Mean'].copy().to_frame()
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_adapter_content(sample, file):
    """Parse fastqc Adapter Content.
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Adapter Content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_per_base_seq_content(sample, file):
    """Parse fastqc Per Base Seq Content
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Per base sequence content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_sequence_length(sample, file):
    """Parse fastqc Sequence Length
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Sequence Length Distribution')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_overrepresented_seq(sample, file):
    """Parse fastqc Overrepresente Sequence
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Overrepresented sequences')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_basic_stats(sample, file):
    """Parse fastqc Basic Stats
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Basic Statistics')
    df = fqc.df.T
    df['sample'] = sample
    return df.set_index('sample')


def parse_fastqc_kmer_content(sample, file):
    """Parse fastqc Kmer Content
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Kmer Content')
    df = fqc.df
    df.reset_index(inplace=True)
    df.set_index('Max Obs/Exp Position', inplace=True)
    splitRanges = split_ranges(df)
    splitRanges.index.name = 'Max Obs/Exp Position'
    splitRanges.reset_index(inplace=True)
    splitRanges['sample'] = sample
    return splitRanges.sort_values(['Sequence', 'Max Obs/Exp Position']).set_index(['sample', 'Sequence'])


def parse_fastqc_per_base_n_content(sample, file):
    """Parse fastqc Per Base N Content
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Per base N content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_per_seq_gc_content(sample, file):
    """Parse fastqc Per Seq GC Content
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Per sequence GC content')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_seq_dup_level(sample, file):
    """Parse fastqc Seq Duplication Level
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    """
    fqc = parse_fastqc(sample, file, field='Sequence Duplication Levels')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()
