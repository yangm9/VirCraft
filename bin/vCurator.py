#!/usr/bin/env python3
import sys
import pandas as pd
from db import suspGenes
#https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=4

def curateV(VirSort2_f,CheckV_f,annot_f):
    '''
    Merge and Filter the contigs according to some certain criteria.
    1) Merge "vs2-pass1/final-viral-score.tsv" and "checkv/contamination.tsv".
    2) Filter the contigs by the empirical screening criteria as follows:
    Keep1: viral_gene >0
    Keep2: viral_gene =0 AND (host_gene =0 OR score >=0.95 OR hallmark >2)
    Manual check: (NOT in Keep1 OR Keep2) AND viral_gene =0 AND host_gene =1 AND length >=10kb
    Discard: the rest
    '''
    df1=pd.DataFrame(pd.read_csv(VirSort2_f,header=0,sep='\t'))
    df1.rename(columns={'seqname':'contig_id'},inplace=1)
    df2=pd.DataFrame(pd.read_csv(CheckV_f,header=0,sep='\t'))
    df=pd.merge(df1,df2,how='left',on='contig_id')
    FiltDF=df[(df['viral_genes']>0)]
    Keep2DF=df[(df['viral_genes']==0)&((df['host_genes']==0)|(df['max_score']>=0.95)|(df['hallmark']>2))]
    Keep2DF,ManuDF=markBySuspGene(Keep2DF,annot_f)
    FiltDF=FiltDF.append(Keep2DF)
    RestDF=df[~df.index.isin(FiltDF.index.tolist())]
    FiltDF=FiltDF.append(RestDF[(RestDF['viral_genes']==0)&(RestDF['host_genes']==1)&(RestDF['length']>10000)])
    return FiltDF,ManuDF

def markBySuspGene(keep2_df,annot_f):
    '''
    Mark the contigs with suspicious gene in dramv-annotate/annotations.tsv for the subsequent manual/semi-autometic curation.
    '''
    df=pd.DataFrame(pd.read_csv(annot_f,header=0,sep='\t'))
    df.rename(columns={'Unnamed: 0':'contig_id'},inplace=1)
    ContigList=[]
    for gene in suspGenes.SuspGeneList:
        SuspDF=df[df['pfam_hits'].str.contains(gene,na=False)]
        tmpList=SuspDF['contig_id'].apply(lambda x:x.split('-')[0].replace('__','||')).tolist()
        ContigList.extend(tmpList)
    ContigList=list(set(ContigList))
    print(ContigList)
    ManuDF=keep2_df[keep2_df['contig_id'].isin(ContigList)]
    keep2_df=keep2_df[~keep2_df['contig_id'].isin(ContigList)]
    return keep2_df,ManuDF

if __name__=='__main__':
    wkdir=sys.argv[1]
    VirSort2_f=wkdir+'/vs2-pass1/final-viral-score.tsv'
    CheckV_f=wkdir+'/checkv/contamination.tsv'
    annot_f=wkdir+'/dramv-annotate/annotations.tsv'
    df,ManuDF=curateV(VirSort2_f,CheckV_f,annot_f)
    df.to_csv(wkdir+'/curation/curated_contigs.xls',index=False,sep='\t')
    ManuDF.to_csv(wkdir+'/curation/manu_curate_contigs.xls',index=False,sep='\t')
