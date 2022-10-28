#!/uer/bin/env python3
import re
import pandas as pd
df=pd.read_csv('all_merged_anno.tpm',sep='\t')
df=df[['Contig', 'V1', 'V2', 'V3', 'V4', 'VR1', 'VR2', 'VR3', 'CB148MP1_CB','CB148MP2_CB', 'CB200MS1_CB', 'CB237MS1_CB', 'CB247MS1_CB','CB267MS1_CB', 'CB267MS2_CB', 'CB900MS1_CB', 'BH0560M_AH', 'BH0960M_AH','BH0995M_AH', 'BH05106M_AH','Order','Family']]
df['Order']=df['Order'].astype('str')
df['Family']=df['Family'].astype('str')
df['Source']=df['Contig'].apply(lambda x:x.split('_')[0])
df['Length']=df['Contig'].apply(lambda x:x.split('_')[4])
df['Contig']=df['Contig'].apply(lambda x:re.sub(r'_length_\d+_cov_\d+\.\d*_\d*','',x))
df['Contig']=df['Contig']+':'+df['Order']+';'+df['Family']
df=df[['Contig', 'V1', 'V2', 'V3', 'CB148MP1_CB','CB148MP2_CB', 'CB200MS1_CB', 'CB237MS1_CB', 'CB247MS1_CB','CB267MS1_CB', 'CB267MS2_CB', 'CB900MS1_CB', 'BH0560M_AH', 'BH0960M_AH','BH0995M_AH', 'BH05106M_AH','Source','Length']]
df.to_csv('all_merged_anno_modi.tpm.xls',sep='\t',index=False)
