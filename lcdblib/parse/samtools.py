import re
import pandas as pd
from collections import OrderedDict

def parse_samtools_stats(sample, file):
    """Parse samtools stats.

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
            if l.startswith('SN'):
                fqs = re.search(r"^SN\s+(.+?):\s+([\d\.]+)\s.*$", l)
                if fqs:
                    if '.' in fqs.group(2):
                        parsed[fqs.group(1)] = float(fqs.group(2))
                    else:
                        parsed[fqs.group(1)] = int(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])
