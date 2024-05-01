#!/usr/bin/env python3

import os
import sys
import pandas as pd

df=pd.read_csv(sys.argv[1])
df['Source']=df['Genome'].apply(lambda x:x.split('_')[0] if not '~' in x else 'ProkaryoticViralRefSeq')
srcs=df['Source'].unique().tolist()
src_num=len(srcs)
bin_dir=os.path.dirname(os.path.realpath(__file__))
cmd_txt='''source "/backup/software/miniconda3/etc/profile.d/conda.sh"
conda activate && conda activate VC-General
Rscript {bin_dir}/venn{}.R'''.format(src_num)
for src in srcs:
    tmp_df=df[df['Source']==src]['VC'].dropna()
    tmp_df.to_csv(src+'.list',sep='\t',index=False,header=False)
    cmd_txt+=' '+src+'.list'
cmd_txt+=' vContact2_venn_diagram.pdf'
os.system(cmd_txt)
