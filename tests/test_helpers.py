import pytest
import pandas as pd
from lcdblib.snakemake import helpers

@pytest.fixture(scope='session')
def patterns():
    patterns = {
        'fastq':   'samples/{sample}/{sample}_R{N}.fastq.gz',
        'cutadapt': 'samples/{sample}/{sample}_R{N}.cutadapt.fastq.gz',
        'bam':     'samples/{sample}/{sample}.cutadapt.bam',
    }

    return patterns

def test_fill_patterns_dict(patterns):
        fill = dict(sample=['U1', 'U2', 'T1', 'T2'], N=[1, 2])
        res = helpers.fill_patterns(patterns, fill)['fastq']

        expected = [
                    'samples/T1/T1_R1.fastq.gz',
                    'samples/T1/T1_R2.fastq.gz',
                    'samples/T2/T2_R1.fastq.gz',
                    'samples/T2/T2_R2.fastq.gz',
                    'samples/U1/U1_R1.fastq.gz',
                    'samples/U1/U1_R2.fastq.gz',
                    'samples/U2/U2_R1.fastq.gz',
                    'samples/U2/U2_R2.fastq.gz',
                   ]

        assert expected == sorted(res)

def test_fill_patterns_dict_zip(patterns):
        fill = dict(sample=['U1', 'T2'], N=[1, 2])
        res = helpers.fill_patterns(patterns, fill, zip)['fastq']

        expected = [
                    'samples/T2/T2_R2.fastq.gz',
                    'samples/U1/U1_R1.fastq.gz',
                   ]

        assert expected == sorted(res)

def test_fill_patterns_df(patterns):
    fill = pd.DataFrame({
        'sample': ['U1', 'U2', 'T1', 'T2'],
        'N': [1, 2, 3, 4],
        })

    res = helpers.fill_patterns(patterns, fill)['fastq']

    expected = [
                'samples/T1/T1_R3.fastq.gz',
                'samples/T2/T2_R4.fastq.gz',
                'samples/U1/U1_R1.fastq.gz',
                'samples/U2/U2_R2.fastq.gz',
               ]

    assert expected == sorted(res)
