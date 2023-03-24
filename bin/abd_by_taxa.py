#!/usr/bin/env python3
#Old Name: specTab4barplot.py

import sys
import pandas as pd

def sumAbundByTaxa(taxa_tpm: str, outdir: str):
    '''
    Sum the abundance by Taxa for each sample.
    input format:
    Contig\tSample1 Abundance\t...\tSampleN Abundance\n
    '''
    df=pd.read_csv(taxa_tpm,sep='\t',header=0)
    samp_num=df.columns.size
    df=df.iloc[:,0:samp_num]
    df['Contig']=df['Contig'].apply(lambda x:x.split(':')[1])
    df['Contig'].replace('nan;nan','Unassigned',inplace=True)
    TaxaList=df['Contig'].unique().tolist()
    df=df.set_index('Contig')
    sumDF=pd.DataFrame()
    for i in range(len(TaxaList)):
        series=df.loc[TaxaList[i]]
        if isinstance(series,pd.DataFrame):
            series=pd.Series(series.sum(),name=TaxaList[i])
        sumDF=pd.concat([sumDF,series.to_frame().T])
        #sumDF=sumDF.append(series) #append will be instead by concat
    sumDF.reset_index(inplace=True)
    sumDF.rename(columns={'index':'Contig'},inplace=True)
    sumDF.fillna(0,inplace=True)
    sumed_tpm=f'{outdir}/all_taxa_sum_abd.xls'
    sumDF.to_csv(sumed_tpm,sep='\t',index=False)
    return 0

if __name__=='__main__':
    if len(sys.argv)==3:
        sumAbundByTaxa(sys.argv[1],sys.argv[2])
    else:
        print(f'Usage: python {sys.argv[0]} <merged_anno_tpm.xls> <outdir>')
