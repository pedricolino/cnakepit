import pandas as pd
from pathlib import Path
import os

def read_sample(data_dir):
    pathlist = Path(data_dir).glob('*1.fastq.gz')
    dic = {}

    for path in pathlist:
        path_str = str(path)
        sample = os.path.basename(path).split('_')[0]
        path_2 =  path_str.replace("1.fastq.gz", "2.fastq.gz")
        if sample not in dic:
            dic[sample] = (sample, path_str, path_2)

    samples = pd.DataFrame.from_dict(dic, orient='index', columns=["samplename", "fq1", "fq2"])
    return samples

if __name__=="__main__":
    print(read_sample("resources/data"))

