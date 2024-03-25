#!/usr/bin/env python3
import sys
import pandas as pd
from collections import Counter

def mostFreqRate(taxa_l):
    countDict=dict(Counter(taxa_l))
    max_k,max_v=max(countDict.items(), key=lambda x:x[1])
    max_rate=max_v/len(taxa_l)+0.0)
    most=''
    if max_rate >= 0.5:
        most=max_k
    else:
        most='NA'
    return most

def classify(blast_sp_f):
    df=pd.read_csv(blast_sp_f,sep='\t')
    df['Contig']=df['QueryID'].apply(lambda x:re.sub(r'_\d+$','',x))
    df.drop(labels=['Taxid','Taxonomy','Levels'],axis=1)
    contigs=df['Contig'].tolist()
    for ctg in contigs:
        MatchedFamilies=df[df['Contig']==ctg]['Family'].tolist()
        family=mostFreqRate(MatchedFamilies)
        print('{ctg}\t\t{family}\n')
    return 0
