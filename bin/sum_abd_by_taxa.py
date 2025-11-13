#!/usr/bin/env python3
#Old Name: specTab4barplot.py
#yangm@idsse.ac.cn

import sys
import pandas as pd

# Input format(taxa_abd):
# taxa_abd file Columns: Contig,Sample_1,Sample_2...Sample_N,Superrealm,...,Species stored abundance information
# Contig\tSample1 Abundance\t...\tSampleN Abundance\n

def idxAbd(taxa_abd: str, taxa: str):  
    taxa_levels = ['Superrealm', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df = pd.read_csv(taxa_abd, sep='\t', header=0)
    df['Class'] = df['Class'].astype('str')
    df['Order'] = df['Order'].astype('str')
    df['Family'] = df['Family'].astype('str')
    if taxa == 'taxa':
        df['Contig'] = 'p__' + df['Phylum'] + ';c__' + df['Class'] + ';o__' +df['Order'] + ';f__' + df['Family']
    elif taxa == 'Class':
        df['Contig'] = 'p__' + df['Phylum'] + ';c__' + df['Class']
    elif taxa == 'Order':
        df['Contig'] = 'c__' + df['Class'] + ';o__' +df['Order']
    elif taxa == 'Family':
        df['Contig'] = 'o__' +df['Order'] + ';f__' + df['Family']
    
    df['Contig'] = df['Contig'].str.replace('nan', '', regex=False)
    df.drop(columns=taxa_levels,inplace=True)
    return df

#Sum the abundance by Taxa for each sample.
def sumAbds(df: str, taxa: str, outdir: str):
    TaxaList = df['Contig'].unique().tolist()
    df = df.set_index('Contig')
    sumDF = pd.DataFrame()
    for i in range(len(TaxaList)):
        series = df.loc[TaxaList[i]]
        if isinstance(series, pd.DataFrame):
            series = pd.Series(series.sum(), name=TaxaList[i])
        sumDF =pd.concat([sumDF, series.to_frame().T])
        #sumDF = sumDF.append(series) #append will be instead by concat
    sumDF.reset_index(inplace=True)
    sumDF.rename(columns={'index': 'Contig'}, inplace=True)
    sumDF.fillna(0, inplace=True)
    sumed_abd = f'{outdir}/summed_v{taxa}_abundance.tsv'
    sumDF.to_csv(sumed_abd, sep='\t', index=False)
    return 0

def sumAbdByTaxa(taxa_abd: str, outdir: str):
    df = idxAbd(taxa_abd, 'taxa')
    sumAbds(df, 'taxa', outdir)
    df = idxAbd(taxa_abd, 'Class')
    sumAbds(df, 'Class', outdir)
    df = idxAbd(taxa_abd, 'Order')
    sumAbds(df, 'Order', outdir)
    df = idxAbd(taxa_abd, 'Family')
    sumAbds(df, 'Family', outdir)
    return 0

if __name__ == '__main__':
    if len(sys.argv) == 3:
        sumAbdByTaxa(sys.argv[1], sys.argv[2])
    else:
        print(f'Usage: python {sys.argv[0]} <merged_anno_abd.tsv> <outdir>')
