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
    if seq_type == 'Gene':
        abd_file_dict = {i: abd_d + '/' + i + '_gene_quant/quant.sf' for i in samples}
    elif seq_type == 'Contig' or 'metabat':
        abd_file_dict = {i: abd_d + '/' + i + '.cov' for i in samples}
    else:
        raise ValueError('Sequence type must be "Contig", "Gene" or "metabat"')
    merged_df = pd.DataFrame()
    seq_dict = {'Contig': 'Contig', 'metabat': 'Contig', 'Gene': 'Gene'}
    for sample_name, file_path in abd_file_dict.items():
        df = pd.read_csv(file_path, sep='\t', header=0)
        if seq_type == 'Gene':
            df = df.iloc[:, [0, 3]]
            df.columns = [seq_dict[seq_type], sample_name]
        elif seq_type == 'metabat':
            df = df.iloc[:, [0, 3, 4]]
            df.columns = [seq_dict[seq_type], sample_name, f'{sample_name}.var']
        else:
            df.columns = [seq_dict[seq_type], sample_name]
        merged_df = df if merged_df.empty else pd.merge(merged_df, df, on=seq_dict[seq_type], how='outer')
    merged_df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) == 5:
        mergeAbundanceFiles(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print(f'Usage: {sys.argv[0]} <contig/gene_abundance_dir> <sample_info.tsv> <output_tsv> <Contig/Gene>')
