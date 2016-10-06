import unittest, pytest
import os
import sys
import re
import yaml
from jsonschema import validate, ValidationError

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

    def test_validate_config(self):
        # TODO: Need to get schema checking and default config setup
        pass

    @unittest.expectedFailure
    def test_wrapper_for(self):
        wrapper_for = helpers.build_wrapper_for(HERE, '../../wrappers')
        result = wrapper_for('cutadapt')
        assert os.path.exists(os.path.join(result, 'wrapper.py'))

    def test_build_params_for(self):
        params_for = helpers.build_params_for(self.config)
        self.assertEqual(params_for('cutadapt', 'extra'), '-a file:adapters.fa --quality-cutoff 20 --minimum-length 25 --overlap 10')

    def test_build_threads_for(self):
        threads_for = helpers.build_threads_for(self.config)
        self.assertEqual(threads_for('cutadapt'), 10)

    @unittest.expectedFailure
    def test_workflow_helper_functions(self):
        wrapper_for, params_for, threads_for = helpers.workflow_helper_functions(self.config, HERE, '../../wrappers')

        # Wrapper For
        result = wrapper_for('cutadapt')
        assert os.path.exists(os.path.join(result, 'wrapper.py'))

        # Params For
        params_for = helpers.build_params_for(self.config)
        self.assertEqual(params_for('cutadapt', 'extra'), '-a file:adapters.fa --quality-cutoff 20 --minimum-length 25 --overlap 10')

        # Threads For
        threads_for = helpers.build_threads_for(self.config)
        self.assertEqual(threads_for('cutadapt'), 10)

    def test_load_sampletable(self):
        from pandas.util.testing import assert_frame_equal
        df = helpers.load_sampletable(self.tempFile)
        assert_frame_equal(df, self.df)

