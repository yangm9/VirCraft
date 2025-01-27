#!/usr/bin/env python3
#Old Name: specTab4barplot.py
#yangm@idsse.ac.cn

import sys
import pandas as pd

#input format:
#Contig\tSample1 Abundance\t...\tSampleN Abundance\n
#Contig format: "contig_id_12345:Caudovirales;Myoviridae"
#taxa:
#taxonomic level, including Order, Family and taxa.
taxa_loc_dict={'Order':0,'Family':1}
def idxAbd(taxa_abd:str,taxa:str,outdir:str):
    df=pd.read_csv(taxa_abd,sep='\t',header=0)
    samp_num=df.columns.size
    df=df.iloc[:,0:samp_num]
    df['Contig']=df['Contig'].apply(lambda x:x.split(':')[1])
    if taxa=='taxa':
        df['Contig'].replace('nan;nan','Unassigned',inplace=True)
    else:
        df['Contig']=df['Contig'].apply(lambda x:x.split(';')[taxa_loc_dict[taxa]])
        df['Contig'].replace('nan','Unassigned',inplace=True)
    return df

#Sum the abundance by Taxa for each sample.
def sumAbds(df:str,taxa:str,outdir:str):
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
    sumed_abd = f'{outdir}/all_{taxa}_sum_abd.tsv'
    sumDF.to_csv(sumed_abd,sep='\t',index=False)
    return 0

def sumAbdByTaxa(taxa_abd:str,outdir:str):
    df=idxAbd(taxa_abd,'taxa',outdir)
    sumAbds(df,'taxa',outdir)
    df=idxAbd(taxa_abd,'Order',outdir)
    sumAbds(df,'Order',outdir)
    df=idxAbd(taxa_abd,'Family',outdir)
    sumAbds(df,'Family',outdir)
    return 0

if __name__=='__main__':
    if len(sys.argv)==3:
        sumAbdByTaxa(sys.argv[1],sys.argv[2])
    else:
        print(f'Usage: python {sys.argv[0]} <merged_anno_abd.tsv> <outdir>')
