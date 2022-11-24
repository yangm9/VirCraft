#!/usr/bin/env python3
import re
import sys
import pandas as pd
import linkTab

if not len(sys.argv)==3:
    

df=pd.read_csv(sys.argv[1],sep='\t')
samp_num=df.columns.size-4
df['Total_Abundance']=df.iloc[:,1:samp_num].sum(axis=1)
df.rename(columns={'Length':'Contig_Length'},inplace=True)
df=df[['Contig','Total_Abundance','Contig_Length','Source']]
pattern=re.compile('\:\w+;\w+$')
df['Contig']=df['Contig'].apply(lambda x:re.sub(pattern,'',x))
sum_tpm_xls=f'{sys.argv[3]}/contig_sum_tpm.xls'
df.to_csv(sum_tpm_xls,sep='\t',index=False)
df=pd.read_csv(sys.argv[2],sep='\t')
df.rename(columns={'contig_id':'Contig'},inplace=True)
pattern=re.compile('_length_\d+_cov_\d+\.\d*_\d*')
df['Contig']=df['Contig'].apply(lambda x:re.sub(pattern,'',x))
df=df[['Contig','checkv_quality']]
qual_id_xls=f'{sys.argv[3]}/contig_quality_summary.xls'
df.to_csv(qual_id_xls,sep='\t',index=False)
len_sum_tpm_qual_xls=f'{sys.argv[3]}/Contig_len_sum_tpm_quality.xls'
linkTab.merge(sum_tpm_xls,qual_id_xls,'left','Contig',len_sum_tpm_qual_xls)
