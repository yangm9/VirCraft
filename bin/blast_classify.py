#!/usr/bin/env python3
import sys
import pandas as pd

def classify(blast_sp_f):
    df=pd.read_csv(blast_sp_f,sep='\t')
    df['Contig']=df['QueryID'].apply(lambda x:re.sub(r'_\d+$','',x))

