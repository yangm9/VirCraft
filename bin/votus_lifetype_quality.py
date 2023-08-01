#!/usr/bin/env python3
import sys
import pandas as pd

def merge_vbType_ckvQual(ckv_qual_f,vb_genome_qual_f):
    vb_genome_qual_df=pd.read_csv(vb_genome_qual_f,sep='\t',low_memory=False)
    vb_genome_qual_df=vb_genome_qual_df.drop('Quality',axis=1).drop_duplicates()
    vb_genome_qual_df.rename(columns={'scaffold':'contig_id'},inplace=True)
    ckv_qual_df=pd.read_csv(ckv_qual_f,sep='\t',low_memory=False)
    df=pd.merge(ckv_qual_df,vb_genome_qual_df,how='left',on='contig_id')
    df['type']=df['type'].fillna('Undetermined')
    return df

if __name__=='__main__':
    df=merge_vbType_ckvQual(sys.argv[1],sys.argv[2])
    df.to_csv(sys.argv[3],sep='\t',index=0)
