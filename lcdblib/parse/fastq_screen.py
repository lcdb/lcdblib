import re
import pandas as pd
from collections import OrderedDict

def parse_fqscreen(sample, file):
    """Parser fastq screen summary table.

    Adapted from multiqc.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.

    """
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(\S+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
            if fqs:
                org = fqs.group(1)
                parsed[(org, 'reads_processed', 'count')] = int(fqs.group(2))
                parsed[(org, 'unmapped', 'count')] = int(fqs.group(3))
                parsed[(org, 'unmapped', 'percent')] = float(fqs.group(4))
                parsed[(org, 'one_hit_one_library', 'count')] = int(fqs.group(5))
                parsed[(org, 'one_hit_one_library', 'percent')] = float(fqs.group(6))
                parsed[(org, 'multiple_hits_one_library', 'count')] = int(fqs.group(7))
                parsed[(org, 'multiple_hits_one_library', 'percent')] = float(fqs.group(8))
                parsed[(org, 'one_hit_multiple_libraries', 'count')] = int(fqs.group(9))
                parsed[(org, 'one_hit_multiple_libraries', 'percent')] = float(fqs.group(10))
                parsed[(org, 'multiple_hits_multiple_libraries', 'count')] = int(fqs.group(11))
                parsed[(org, 'multiple_hits_multiple_libraries', 'percent')] = float(fqs.group(12))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])
