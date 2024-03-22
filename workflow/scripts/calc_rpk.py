import pandas as pd
import sys
from pathlib import Path

stat_file = snakemake.input[0]
df = pd.read_csv(stat_file, sep='\t')
rpk = df.iloc[:, 2] / df.iloc[:, 1] * 1000
df.insert(4, None, rpk)
df.to_csv(str(Path(stat_file)) + "_aug", sep='\t', index=False, header=None)
