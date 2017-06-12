from io import StringIO
import pandas as pd

def parse_picardCollect_summary(sample, file):
    """Parser for picard collectRNAMetrics summary.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    with open(file, 'r') as fh:
        for l in fh:
            if l.startswith('#'):
                continue
            if l.startswith('PF_BASES'):
                parsed = l
                parsed += next(fh)
                break

        if len(parsed) == 0:
            return None
        else:
            df = pd.read_csv(StringIO(parsed), sep='\t')
            df.index = [sample]
            return df


def parse_picardCollect_hist(sample, file):
    """Parser for picard collectRNAMetrics summary.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    parsed = ''
    with open(file, 'r') as fh:
        for l in fh:
            if l.startswith('#'):
                continue
            if l.startswith('normalized'):
                parsed = l
                while True:
                    try:
                        parsed += next(fh)
                    except StopIteration:
                        break
                break

        if len(parsed) == 0:
            return None
        else:
            df = pd.read_csv(StringIO(parsed), sep='\t', index_col=0).T
            df.index = [sample]
            return df


def parse_picard_markduplicate_metrics(sample, file):
    """Parser for picard markduplicates.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    df = pd.read_csv(file, sep='\t', comment='#')
    df['sample'] = sample
    return df.set_index('sample')
