#!/usr/bin/env python
import unittest
import pandas as pd
from pandas.util.testing import assert_frame_equal

from lcdblib.parse import fastqc

BLOCK = """>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	33.325882537780295	34.0	34.0	34.0	31.0	34.0
2	33.44737950739895	34.0	34.0	34.0	31.0	34.0
3	33.446659779075645	34.0	34.0	34.0	31.0	34.0
4	36.67489445866323	37.0	37.0	37.0	37.0	37.0
5	36.658260994356844	37.0	37.0	37.0	35.0	37.0
6	36.65553688182567	37.0	37.0	37.0	35.0	37.0
7	36.65730974503444	37.0	37.0	37.0	35.0	37.0
8	36.62655009906297	37.0	37.0	37.0	35.0	37.0
9	38.53114508354548	39.0	39.0	39.0	38.0	39.0
10	38.47594460091534	39.0	39.0	39.0	37.0	39.0
"""

BLOCK_DF = pd.DataFrame({'Base': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                         'Mean': [33.325882537780295, 33.44737950739895, 33.446659779075645,
                                  36.67489445866323, 36.658260994356844, 36.65553688182567,
                                  36.65730974503444, 36.62655009906297, 38.53114508354548,
                                  38.47594460091534],
                         'Median': [34.0, 34.0, 34.0, 37.0, 37.0, 37.0, 37.0, 37.0, 39.0, 39.0],
                         'Lower Quartile': [34.0, 34.0, 34.0, 37.0, 37.0, 37.0, 37.0, 37.0, 39.0, 39.0],
                         'Upper Quartile': [34.0, 34.0, 34.0, 37.0, 37.0, 37.0, 37.0, 37.0, 39.0, 39.0],
                         '10th Percentile': [31.0, 31.0, 31.0, 37.0, 35.0, 35.0, 35.0, 35.0, 38.0, 37.0],
                         '90th Percentile': [34.0, 34.0, 34.0, 37.0, 37.0, 37.0, 37.0, 37.0, 39.0, 39.0]})

BLOCK_DF.set_index('Base', inplace=True)


class TestFastqc(unittest.TestCase):
    def setUp(self):
        self.filename = 'tests/data/fastqc_data.txt'
        self.id = 'test'

    def tearDown(self):
        pass

    def test_FastQCBlock(self):
        fq = fastqc.FastQCBlock(BLOCK.split('\n'))
        assert_frame_equal(fq.df, BLOCK_DF[fq.df.columns])
        
    def test_fileParse(self):
        fq = fastqc.FastQC('test', 'tests/data/fastqc_data.txt')
        assert_frame_equal(fq['Per base sequence quality'].df, 
                           BLOCK_DF[fq['Per base sequence quality'].df.columns])




if __name__ == '__main__':
    unittest.main()

