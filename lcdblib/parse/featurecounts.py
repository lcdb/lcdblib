import re
import pandas as pd
from collections import OrderedDict

def parse_featureCounts_counts(sample, file):
    """Parser for featurecounts counts table.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    df = pd.read_csv(file, sep='\t', comment='#')
    df.columns = ['FBgn', 'chr', 'start', 'end', 'strand', 'length', 'count']
    df['sample'] = sample
    df.set_index(['sample', 'FBgn'], inplace=True)
    return df['count']


def parse_featureCounts_summary(sample, file):
    """Parse featurecounts summary table

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?)\s+(\d+)$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])
