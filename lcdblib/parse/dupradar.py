import pandas as pd

def parse_dupradar(sample, file):
    """Parser for dupradar.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    df = pd.read_csv(file, sep='\t', index_col=0)
    df.columns = ['FBgn', 'geneLength', 'allCountsMulti', 'filteredCountsMulti', 'dupRateMulti',
                  'dupsPerIdMulti', 'RPKMulti', 'RPKMMulti', 'allCounts', 'filteredCounts',
                  'dupRate', 'dupsPerId', 'RPK', 'RPKM', 'mhRate']
    df['sample'] = sample
    return df.set_index(['sample', 'FBgn'])
