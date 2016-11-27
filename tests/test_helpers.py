import unittest, pytest
import os
import sys
import re
import yaml

import pandas as pd

from lcdblib.snakemake import helpers

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(HERE, '../'))


class TestHelperFunctions(unittest.TestCase):

    def setUp(self):

        # Example config
        self.config = {'sampletable': 'pasilla_sampletable.tsv',
                       'sample_dir': 'pasilla',
                       'assembly': 'dm6',
                       'data_dir': '/data/LCDB/references-test',
                       'shell_prefix': "source lcdb-workflows.bashrc; ",
                       'rules': {
                           'cutadapt': {
                               'params': {
                                   'extra': '-a file:adapters.fa --quality-cutoff 20 --minimum-length 25 --overlap 10'
                               },
                               'threads': 10
                           }
                       }
                       }

        # Build basic sample table
        from tempfile import mktemp
        self.tempFile = mktemp()

        self.df = pd.DataFrame({'sampleID': ['U1', 'U2', 'T1', 'T2'], 'sex': ['M', 'M', 'M', 'M'], 'rep': [1, 2, 1, 2], })
        self.df.set_index('sampleID', inplace=True)
        self.df.to_csv(self.tempFile, sep='\t')

    def tearDown(self):
        os.remove(self.tempFile)

