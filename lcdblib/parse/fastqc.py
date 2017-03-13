#!/usr/bin/env python
""" Quick and dirty Fastqc parser """
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
