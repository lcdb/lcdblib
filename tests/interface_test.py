import unittest
import sys
import os
import yaml

# This may not be needed if tests are always ran from project folder, but
# guarantees that paths are correct.
HERE = os.path.dirname(os.path.realpath(__file__))
os.chdir(HERE)

configYaml = """
sampletable: test_sampletable.tsv
sample_dir: pasilla
assembly: dm6
data_dir: references

rawLevel: 'pasilla/{sampleID}/{sampleID}_{treatment}_{replicate}_0001_R1'
runLevel: 'pasilla_sample/{sampleID}/{sampleID}_{treatment}_{replicate}_R1'
sampleLevel: 'pasilla_sample/{sampleID}/{sampleID}_{treatment}'
aggLevel: 'pasilla_agg/{treatment}'
aggLevel2:
    - '{treatment}'
"""

test_sampletable = """\
sampleID	treatment	replicate
treated1	treated	1
treated2	treated	2
untreated1	untreated	1
untreated2	untreated	2
"""
sampletable = os.path.join(os.path.dirname(__file__), 'test_sampletable.tsv')

samples = [
    {'sampleID': 'treated1', 'replicate': '1', 'treatment': 'treated'},
    {'sampleID': 'treated2', 'replicate': '2', 'treatment': 'treated'},
    {'sampleID': 'untreated1', 'replicate': '1', 'treatment': 'untreated'},
    {'sampleID': 'untreated2', 'replicate': '2', 'treatment': 'untreated'}
]

patterns = [
        '{runLevel}.fastq.gz',
        '{runLevel}.sort.bam',
        '{sampleLevel}.merged.sort.bam',
        'otherfile.txt'
        ]


class TestSampleHandler(unittest.TestCase):
    def setUp(self):
        from lcdblib.snakemake.interface import SampleHandler

        with open(sampletable, 'w') as fout:
            fout.write(test_sampletable)

        # Load config
        self.config = yaml.load(configYaml)

        # Initialize the class
        self.SH = SampleHandler(self.config)

    def tearDown(self):
        os.unlink('test_sampletable.tsv')

    def test_init(self):
        self.assertEqual(self.config, self.SH.config)
        self.assertEqual(samples, self.SH.samples)

    def test_find_level_raw(self):
        prefix = self.config['rawLevel'].format_map(samples[0])
        self.assertEqual('rawLevel', self.SH.find_level(prefix)[0])

    def test_find_level_run(self):
        prefix = self.config['runLevel'].format_map(samples[0])
        self.assertEqual('runLevel', self.SH.find_level(prefix)[0])

    def test_find_level_sample(self):
        prefix = self.config['sampleLevel'].format_map(samples[0])
        self.assertEqual('sampleLevel', self.SH.find_level(prefix)[0])

    def test_find_level_bad_pattern(self):
        prefix = '../../test/{sampleID}/{sampleID}/{sampleID}'.format_map(samples[0])
        with self.assertRaises(ValueError):
            self.SH.find_level(prefix)

    def test_find_sample_raw(self):
        prefix = self.config['rawLevel'].format_map(samples[0])
        _samples = self.SH.find_sample(self.SH.raw, prefix)
        self.assertEqual(_samples[0], samples[0])

    def test_find_sample_run(self):
        prefix = self.config['runLevel'].format_map(samples[0])
        _samples = self.SH.find_sample(self.SH.run, prefix)
        self.assertEqual(_samples[0], samples[0])

    def test_find_sample_sample(self):
        prefix = self.config['sampleLevel'].format_map(samples[0])
        _samples = self.SH.find_sample(self.SH.sample, prefix)
        self.assertEqual(_samples[0], samples[0])

    def test_make_input_raw(self):
        wildcards = {'prefix': self.config['rawLevel'].format_map(self.SH.samples[0])}
        _input = self.SH.make_input(suffix='.fastq', agg=False)
        self.assertEqual(_input(wildcards), [wildcards['prefix'] + '.fastq'])

    def test_make_input_run(self):
        wildcards = {'prefix': self.config['runLevel'].format_map(self.SH.samples[0])}
        _input = self.SH.make_input(suffix='.fastq', agg=False)
        self.assertEqual(_input(wildcards), [wildcards['prefix'] + '.fastq'])

    def test_make_input_sample(self):
        wildcards = {'prefix': self.config['sampleLevel'].format_map(self.SH.samples[0])}
        _input = self.SH.make_input(suffix='.fastq', agg=False)
        self.assertEqual(_input(wildcards), [wildcards['prefix'] + '.fastq'])

    def test_make_input_run_agg(self):
        wildcards = {'prefix': self.config['runLevel'].format_map(self.SH.samples[0])}
        rawInput = self.config['rawLevel'].format_map(self.SH.samples[0])
        _input = self.SH.make_input(suffix='.fastq', agg=True)
        self.assertEqual(_input(wildcards), [rawInput + '.fastq'])

    def test_make_input_sample_agg(self):
        wildcards = {'prefix': self.config['sampleLevel'].format_map(self.SH.samples[0])}
        runInput = self.config['runLevel'].format_map(self.SH.samples[0])
        _input = self.SH.make_input(suffix='.fastq', agg=True)
        self.assertEqual(_input(wildcards), [runInput + '.fastq'])

    def test_make_input_sample_agg_additional_prefix(self):
        wildcards = {
                'prefix': self.config['sampleLevel'].format_map(self.SH.samples[0]),
                'midfix': '.cutadapt.bt2',
                'suffix': '.unique.sort.dedup.bam'}

        runInput = self.config['runLevel'].format_map(self.SH.samples[0]) + '.cutadapt.bt2.unique.sort.dedup.bam'

        _input = self.SH.make_input(prefix='prefix', midfix='midfix', suffix='suffix', agg=True)
        self.assertEqual(_input(wildcards), [runInput])

    def test_build_targets(self):
        targets = [
                'pasilla_sample/treated1/treated1_treated_1_R1.fastq.gz',
                'pasilla_sample/treated2/treated2_treated_2_R1.fastq.gz',
                'pasilla_sample/untreated1/untreated1_untreated_1_R1.fastq.gz',
                'pasilla_sample/untreated2/untreated2_untreated_2_R1.fastq.gz',
                'pasilla_sample/treated1/treated1_treated_1_R1.sort.bam',
                'pasilla_sample/treated2/treated2_treated_2_R1.sort.bam',
                'pasilla_sample/untreated1/untreated1_untreated_1_R1.sort.bam',
                'pasilla_sample/untreated2/untreated2_untreated_2_R1.sort.bam',
                'pasilla_sample/treated1/treated1_treated.merged.sort.bam',
                'pasilla_sample/treated2/treated2_treated.merged.sort.bam',
                'pasilla_sample/untreated1/untreated1_untreated.merged.sort.bam',
                'pasilla_sample/untreated2/untreated2_untreated.merged.sort.bam',
                'otherfile.txt'
                ]
        self.assertEqual(set(self.SH.build_targets(patterns)), set(targets))
