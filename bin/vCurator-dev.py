#!/usr/bin/env python3
import sys
import pandas as pd
from db import suspGenes

#Reference DOI: dx.doi.org/10.17504/protocols.io.bwm5pc86

def autoCurate(VirSort2_f,CheckV_f,anno_f):
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
    return ManuList,ManuDF,AnnoDF

#ContigAnnoDF.iloc[0,[8,10,16,23,24,26]]
def calcAnno(ctg_anno_df):
    AnnoDbIdx=[8,10,14,16,23,24,26]
    ContigInfo={
        'anno':0,'ko':0,'viral':0,'viral_gt_100':0,
        'peptidase':0,'pfam':0,'cazy':0,'vogdb':0
    }
    row_num=ctg_anno_df.shape[0]
    for i in range(row_num):
        SubSeries=ctg_anno_df.iloc[i,AnnoDbIdx].fillna(0)
        if SubSeries.any():
            ContigInfo['anno']+=1
            if SubSeries['viral_id']:
                ContigInfo['viral']+=1
                if SubSeries['viral_bitScore'] >= 100:
                    ContigInfo['viral_gt_100']+=1
    return ContigInfo

def contigCurate(contig_info_dict):
    if contig_info_dict['anno']==0: return 'N'
    viral_rate=contig_info_dict['viral']/contig_info_dict['anno']+0.0
    if viral_rate < 0.5:
        return ''
    elif viral_rate >= 0.5:
        viral_gt_100_rate=contig_info_dict['viral_gt_100']/contig_info_dict['viral']+0.0
        if viral_gt_100_rate >= 0.5:
            return 'P'
        else:
            return 'VUS'
    else:
        return ''

def manuCurate(manu_f,anno_f):
    '''
    The classification of viral curation includes Positive (P), Likely Positive(LP), Virus of Uncertain Significance (VUS), Likely Negative (LN) and Negative (N). 
    '''
    ManuList,ManuDF,ManuAnnoDF=filtAnnoForManu(manu_f,anno_f)
    manu_anno_f=wkdir+'/curation/manu_curate_anno.xls'
    ManuAnnoDF.to_csv(manu_anno_f,index=False,sep='\t')
    ManuDF.set_index('contig_id',inplace=True)
    ManuDF['curation']=''
    for contig in ManuList:
        ContigAnnoDF=ManuAnnoDF[ManuAnnoDF['contig_id']==contig]
        ContigInfo=calcAnno(ContigAnnoDF)
        judge=contigCurate(ContigInfo)
        ManuDF.loc[contig,'curation']=judge
    ManuDF.reset_index(inplace=True)
    return ManuDF

if __name__=='__main__':
    wkdir=sys.argv[1]
    VirSort2_f=wkdir+'/vs2-pass1/final-viral-score.tsv'
    CheckV_f=wkdir+'/checkv/contamination.tsv'
    anno_f=wkdir+'/dramv-annotate/annotations.tsv'
    AutoDF,ManuDF=autoCurate(VirSort2_f,CheckV_f,anno_f)
    autu_curated_f=wkdir+'/curation/autu_curated_contigs.xls'
    manu_curate_f=wkdir+'/curation/manu_curate_contigs.xls'
    AutoDF.to_csv(autu_curated_f,index=False,sep='\t')
    ManuDF.to_csv(manu_curate_f,index=False,sep='\t')
    CuratedDF=manuCurate(manu_curate_f,anno_f)
    manu_curated_f=wkdir+'/curation/manu_curated_contigs.xls'
    CuratedDF.to_csv(manu_curated_f,index=False,sep='\t')

'''
    Step 5: Manual curation
    For those in “manual check” category, you can look through their annotations in "dramv-annotate/annotations.tsv", in which each gene of every contig is a line and has annotation from multiple databases. This step is hard to standardize, but below are some criteria based on our experience.
    Criteria for calling a contig viral:
    - Structural genes, hallmark genes, depletion in annotations or enrichment for hypotheticals (~10% genes having non-hypothetical annotations)
    + Lacking hallmarks but >=50% of annotated genes hit to a virus and at least half of those have viral bitcore >100 and the contig is <50kb in length
    - Provirus: Integrase/recombinase/excisionase/repressor, enrichment of viral genes on one side
    - Provirus: “break” in the genome: gap between two genes corresponding to a strand switch, higher coding density, depletion in annotations, and an enrichment for phage genes on one side
    - Few annotations only ~1-3 genes, but with at least half hitting to viruses, and where the genes hitting cells have a bitscore no more than 150% that of the viral bitscores and/or viral bitscores are >100
    - LPS (lipopolysaccharide) looking regions if also has very strong hits to viral genes bitscore > 100
    Criteria for callling a contig non-viral:
    - >3x cellular like genes than viral, nearly all genes annotated, no genes hitting to only viruses and no viral hallmark genes
    - Lacking any viral hallmark genes and >50kb
    - Strings of many obvious cellular genes, with no other viral hallmark genes. Examples encountered in our benchmarking include 1) CRISPR Cas, 2) ABC transporters, 3) Sporulation proteins, 4) Two-component systems, 5) Secretion system. Some of these may be encoded by viruses, but are not indicative of a viral contig without further evidence.
    - Multiple plasmid genes or transposases but no clear genes hitting only to viruses
    - Few annotations, only ~1-3 genes hitting to both viruses and cellular genes but with stronger bitscores for the cellular genes.
    - LPS looking regions if no strong viral hits. Enriched in genes commonly associated with Lipopolysaccharide or LPS, such as epimerases, glycosyl transferases, acyltransferase, short-chain dehydrogenase/reductase, dehydratase
    - Genes annotated as Type IV and/or Type VI secretion system surrounded by non-viral genes
    - Few annotations, only ~1-3 genes all hitting to cellular genes (even if bitscore <100) with no viral hits
    Lastly, user beware that any provirus boundary predicted by VirSorter 2 and/or checkV is an approximate estimate only (calling “ends” is quite a challenging problem in prophage discovery), and needs to be manually inspected carefully too, especially for AMG studies.
'''
