#!/usr/bin/env python3
import re
import sys
import pandas as pd

df=pd.read_csv(sys.argv[1],sep='\t')
df['Order']=df['Order'].astype('str')
df['Family']=df['Family'].astype('str')
df['Source']=df['Contig'].apply(lambda x:x.split('_')[0] if re.search(r'_',x) else x)
#df['Length']=df['Contig'].apply(lambda x:x.split('_')[4])
df['Contig']=df['Contig'].apply(lambda x:re.sub(r'_length_\d+_cov_\d+\.\d*_\d*','',x))
df['Contig']=df['Contig']+':'+df['Order']+';'+df['Family']
df.drop(columns=['Percent_of_votes','Percent_of_votes.1'],inplace=True)
df.to_csv(sys.argv[2],sep='\t',index=False)
