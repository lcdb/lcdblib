import re
from io import StringIO
import pandas as pd
from collections import OrderedDict

def parse_atropos(sample, file):
    """Parse atropos.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.

    """
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        block = None
        subBlock = None
        cnts = {}
        for l in fh:
            l = l.replace(',', '')

            if l.startswith('==='):
                block = re.search(r"^=== (.+) ===$", l).group(1)
            else:
                # The summary block is a little different between SE and PE reads
                if block == 'Summary':
                    fqs = re.search(r"^(\w+[\s\w\(\)-]+?):\s+([\d\.]+)\s.*$", l)
                    fqs2 = re.search(r"^\s+([\w\s-]+?):\s+([\d\.]+)\s.*$", l)

                    if fqs:
                        subBlock = fqs.group(1)
                        parsed[subBlock] = int(fqs.group(2))
                    elif fqs2:
                        read = '_'.join([subBlock, fqs2.group(1)])
                        parsed[read] = int(fqs2.group(2))
                else:
                    fqs = re.search(r"^.*Trimmed:\s+(\d+)\s.*$", l)
                    if fqs:
                        key = 'Number {} trimmed'.format(block)
                        parsed[key] = int(fqs.group(1))
                    elif l.startswith('length'):
                        # This will pull out the length count tables and make dataframes
                        cnts[block] = l
                        try:
                            while True:
                                l = next(fh)
                                if l.startswith('\n'):
                                    break
                                cnts[block] += l
                        except StopIteration:
                            pass
                        cnts[block] = pd.read_table(StringIO(cnts[block]))
                        cnts[block]['adapter'] = block
                        cnts[block]['sample'] = sample
                        cnts[block].set_index(['sample', 'adapter', 'length'], inplace=True)

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[sample])

            if [x for x in df.columns if 'Read 1' in x]:
                # PE
                df['pct_read1_adapters'] = df['Total read pairs processed_Read 1 with adapter'] / df['Total read pairs processed'] * 100
                df['pct_read2_adapters'] = df['Total read pairs processed_Read 2 with adapter'] / df['Total read pairs processed'] * 100
            else:
                # SE
                df['pct_read1_adapters'] = df['Reads with adapters'] / df['Total reads processed'] * 100

            return df, pd.concat(cnts.values())
