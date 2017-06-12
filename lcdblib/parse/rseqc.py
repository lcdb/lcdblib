import pandas as pd
import re
from collections import OrderedDict

def parse_inferExperiment(sample, file):
    """Parse rseqc infer expeirment.
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
            fqs = re.search(r"^(.+?):\s+([\d\.]+)$", l)
            if fqs:
                parsed[fqs.group(1)] = float(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_geneBodyCoverage(sample, file):
    """Parse rseqc genebody coverage.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    with open(file, 'r') as fh:
        lines = fh.readlines()
        header = lines[0].strip().split('\t')[1:]
        values = lines[1].strip().split('\t')[1:]
        parsed = OrderedDict()
        for k, v in zip(header, values):
            parsed[int(k)] = float(v)
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_bamStat(sample, file):
    """Parse rseqc bam stat.

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
            fqs = re.search(r"^(.+?):\s*(\d+)$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_tin(sample, file):
    """Parse rseqc tin.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    with open(file, 'r') as fh:
        lines = fh.readlines()
        header = lines[0].strip().split('\t')[1:]
        values = lines[1].strip().split('\t')[1:]
        parsed = OrderedDict()
        for k, v in zip(header, values):
            parsed[k] = float(v)
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])
