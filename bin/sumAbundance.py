#!/usr/bin/env python3
import re
import sys
import pandas as pd
import linkTab

df=pd.read_csv(sys.argv[1],sep='\t')
samp_num=df.columns.size
df['Total_Abundance']=df.iloc[:,1:samp_num].sum(axis=1)
df=df[['Contig','Total_Abundance']]
pattern=re.compile('\:\w+;\w+$')
df['Contig']=df['Contig'].apply(lambda x:re.sub(pattern,'',x))
df.to_csv(sys.argv[2],sep='\t',index=False)
