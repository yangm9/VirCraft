#!/usr/bin/env python3
import sys
import pandas as pd
from db import suspGenes

#Reference DOI: dx.doi.org/10.17504/protocols.io.bwm5pc86

def curateV(VirSort2_f,CheckV_f,anno_f):
    '''
    Step 4: Screening based on viral and host gene counts, score, hallmark gene counts, and contig length
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
    FiltDF=df[(df['viral_genes']>0)] #Get Keep1
    Keep2DF=df[(df['viral_genes']==0)&((df['host_genes']==0)|(df['max_score']>=0.95)|(df['hallmark']>2))]
    Keep2DF,ManuDF=markBySuspGene(Keep2DF,anno_f) #keep2
    FiltDF=FiltDF.append(Keep2DF) 
    RestDF=df[~df.index.isin(FiltDF.index.tolist())]
    ManuDF=ManuDF.append(RestDF[(RestDF['viral_genes']==0)&(RestDF['host_genes']==1)&(RestDF['length']>10000)])
    return FiltDF,ManuDF

def markBySuspGene(keep2_df,anno_f):
    '''
    Mark the contigs with suspicious gene in dramv-annotate/annotations.tsv for the subsequent manual/semi-autometic curation.
    DRAMv annotation screening
    There are some genes that are common in both viruses and hosts (e.g.  Polyliposaccharides [LPS] related) and mobile element, which can cause false positives in the above "Keep2" category. Thus we want to be cautious with contigs with these genes. We have compiled a list of "suspicious" genes in this link (https://bitbucket.org/MAVERICLab/virsorter2-sop/raw/03b8f28bee979e2b7fd99d7375d915c29c938339/resource/suspicious-gene.list). You can subset the DRAMv table using contigs in the "Keep2" category, and screen for the "suspicious" genes in the subset DRAMv table (ignore case, e.g. use "-i" option for "grep"),  and then put contigs with those genes in the "Manual check" category.
    '''
    df=pd.DataFrame(pd.read_csv(anno_f,header=0,sep='\t'))
    df.rename(columns={'Unnamed: 0':'contig_id'},inplace=1)
    ContigList=[]
    for gene in suspGenes.SuspGeneList:
        SuspDF=df[df['pfam_hits'].str.contains(gene,na=False)]
        tmpList=SuspDF['contig_id'].apply(lambda x:x.split('-')[0].replace('__','||')).tolist()
        ContigList.extend(tmpList)
    ContigList=list(set(ContigList))
    #print(ContigList)
    ManuDF=keep2_df[keep2_df['contig_id'].isin(ContigList)]
    keep2_df=keep2_df[~keep2_df['contig_id'].isin(ContigList)]
    return keep2_df,ManuDF

def filtAnnoForManu(manu_f,anno_f):
    ManuDF=pd.read_csv(manu_f,header=0,sep='\t')
    AnnoDF=pd.read_csv(anno_f,header=0,sep='\t')
    AnnoDF.rename(columns={'Unnamed: 0':'gene_id'},inplace=1)
    AnnoDF['contig_id']=AnnoDF['gene_id'].apply(lambda x:x.split('-')[0].replace('__','||'))
    ManuList=ManuDF['contig_id'].tolist()
    AnnoDF=AnnoDF[AnnoDF['contig_id'].isin(ManuList)]
    return AnnoDF

if __name__=='__main__':
    wkdir=sys.argv[1]
    VirSort2_f=wkdir+'/vs2-pass1/final-viral-score.tsv'
    CheckV_f=wkdir+'/checkv/contamination.tsv'
    anno_f=wkdir+'/dramv-annotate/annotations.tsv'
    AutoDF,ManuDF=curateV(VirSort2_f,CheckV_f,anno_f)
    autu_curated_f=wkdir+'/curation/autu_curated_contigs.xls'
    manu_curate_f=wkdir+'/curation/manu_curate_contigs.xls'
    AutoDF.to_csv(autu_curated_f,index=False,sep='\t')
    ManuDF.to_csv(manu_curate_f,index=False,sep='\t')
    ManuAnnoDF=filtAnnoForManu(manu_curate_f,anno_f)
    manu_anno_f=wkdir+'/curation/manu_curate_anno.xls'
    ManuAnnoDF.to_csv(manu_anno_f,index=False,sep='\t')
