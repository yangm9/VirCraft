#!/usr/bin/env python3
import os
import sys
import pandas as pd

def getSampList(samp_info_tsv: str):
    samples = []
    SAMPLEINFO = open(samp_info_tsv)
    SAMPLEINFO.readline()
    for line in SAMPLEINFO:
        items = line.strip().split('\t')
        sample = items[0]
        samples.append(sample)
    SAMPLEINFO.close()
    return samples

def mergeAbundanceFiles(abd_d, samp_info_tsv, output_tsv, seq_type='Contig'):
    samples = getSampList(samp_info_tsv)
    if seq_type == 'Contig':
        abd_file_dict = {i: abd_d + '/' + i + '.cov' for i in samples}
    elif seq_type == 'Gene':
        abd_file_dict = {i: abd_d + '/' + i + '_gene_quant/quant.sf' for i in samples}
    else:
        raise ValueError("Sequence type must be 'Contig' or 'Gene'")
    merged_df = pd.DataFrame()
    for sample_name, file_path in abd_file_dict.items():
        df = pd.read_csv(file_path, sep='\t', header=0)
        df = df.iloc[:, [0, 3]] if seq_type == 'Gene' else df
        df.columns = [seq_type, f"{sample_name}"]
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=seq_type, how="outer")
    merged_df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) == 5:
        mergeAbundanceFiles(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print(f'Usage: {sys.argv[0]} <contig/gene_abundance_dir> <sample_info.tsv> <output_tsv> <Contig/Gene>')
