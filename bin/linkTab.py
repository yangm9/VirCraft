#!/usr/bin/env python3
import sys
import pandas as pd

def merge(tab1_f,tab2_f,linkType_s,title,out_f):
    df1=pd.DataFrame(pd.read_csv(tab1_f,header=0,sep='\t',low_memory=False))
    df2=pd.DataFrame(pd.read_csv(tab2_f,header=0,sep='\t',low_memory=False))
    df12=pd.merge(df1,df2,how=linkType_s,on=title)
    df12.to_csv(out_f,sep='\t',index=0)
    return 0

def multiple(samp_1,suffix,linkType_s,title,out_f):
    n=0
    name=''
    prefix=''
    for samp in samp_l:
        tab2=samp+suffix
        if n==1:
            prefix=samp
        elif n==2:
            tab1=name+suffix
            tab3=f'{name}{samp}.merged.tmp'
            merge(tab1,tab2,linkType_s,title,tab3)
        else:
            tab1=f'{name}.merged.tmp'
            tab3=''
            merge(tab1,tab2,linkType_s,title,tab3)
        name+=samp
        n+=1
    return 0

if __name__=='__main__':
    if len(sys.argv)==6:
        merge(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    else:
        print('Usage: '+sys.argv[0]+' <table_1> <table_2> <link_type: left/right/inner/outer/> <the upper left corner value> <output_file>')
