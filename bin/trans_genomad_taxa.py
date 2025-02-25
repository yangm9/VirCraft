#!/usr/bin/env python3
import sys
import pandas as pd

def extract_order_family(taxonomy: str):
    if ';' not in taxonomy:
        return 'Unassigned;Unassigned'
    Class,order,family = taxonomy.split(';')[4:7]
    if not family:
        family = 'Unassigned'
    if not order:
        if Class:
            order = f'order_from_{Class}'
        else:
            order = 'Unassigned'
    return f'{order};{family}'

def get_seq_taxa(genomad_f: str, seq_taxa_f: str):
    df = pd.read_csv(genomad_f, sep='\t')
    df_taxa = df[['seq_name', 'taxonomy']]
    df_taxa['OrderFamily'] = df_taxa['taxonomy'].apply(lambda x: extract_order_family(x))
    df_taxa[['Order', 'Family']] = df_taxa['OrderFamily'].str.split(';', expand=True, n=1)
    df_taxa.drop(['taxonomy', 'OrderFamily'], axis=1, inplace=True)
    df_taxa.rename(columns={'seq_name': 'Sequence_ID'}, inplace=True)
    df_taxa.to_csv(seq_taxa_f, sep='\t', index=False)
    return seq_taxa_f

if __name__ == '__main__':
    if len(sys.argv) == 3:
        get_seq_taxa(sys.argv[1], sys.argv[2])
    else:
        print(f'Usage: {sys.argv[0]} <genomad_summary.tsv> <sequence.taxa.txt>')

