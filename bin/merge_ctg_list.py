#!/usr/bin/env python3
import sys
import pandas as pd

PathDict={
    'virsorter2':'vs2-pass1/final-viral-score.tsv',
    'vibrant':'VIBRANT_/VIBRANT_results_scaffolds/VIBRANT_machine_scaffolds.tsv',
    'deepvirfinder':'deepvirfinder'
}

FiltDict={
    'virsorter2':'dsDNAphage>=0.9',
    'vibrant':"prediction=='virus'",
    'deepvirfinder':'score > 0.9 and pvalue<=0.1'
}

NameDict={
    'virsorter2':'seqname',
    'vibrant':'scaffold',
    'deepvirfinder':'name'
}

FastaDict={
    'virsorter2':'vs2-pass1/final-viral-combined.fa',
    'vibrant':'',
    'deepvirfinder':''
}

def vCtgFilt(name,tool,wkdir):
    result=wkdir+'/'+PathDict[tool]
    if tool=='deepvirfinder':
        file_name=os.listdir(result)[0]
        result+='/'+file_name
    elif tool=='vibrant':
        pass
    df=pd.read_csv(result,sep='\t')
    df=df.query(FiltDict[tool])
    return df

def vCtgMerge(wkdir):
    full_ctgs=[]
    vs2_partial_ctgs=[]
    vb_partial_ctgs=[]
    for tool in FiltDict.keys():
        df=vCtgFilt(tool,wkdir)
        if tool=='virsorter2':
            df['seqname']=df['seqname'].apply(lambda x:x.rsplit('_',1)[0] if x.endswith('full') or x.endswith('lt2gene') else x)
            full_tmps=df.query('not seqname.str.contains("partial")')['seqname'].tolist()
            partial_tmps=df.query('seqname.str.contains("partial")')['seqname'].tolist()
            full_ctgs.extend(full_tmps)
            vs2_partial_ctgs.extend(partial_tmps)
        elif tool=='vibrant':
            full_tmps=df.query('not scaffold.str.contains("fragment")')['scaffold']
            partial_tmps=df.query('scaffold.str.contains("fragment")')['scaffold']
            full_ctgs.extend(full_tmps)
            vb_partial_ctgs.extend(partial_tmps)
        else:
            full_ctgs.extend(df[NameDict[tool]].tolist())
    full_ctgs=list(set(full_ctgs))
    partial_ctgs=list(set(partial_ctgs))
    return full_ctgs,vs2_partial_ctgs,vb_partial_ctgs

def listToFile(list_l,list_f):
    LIST=open(list_f,'w')
    for ctg in list_l:
        LIST.write(f'{ctg}\n')
    return 0

def ctgList(wkdir):
    full_ctgs,vs2_partial_ctgs,vb_partial_ctgs=vCtgMerge(wkdir)
    full_ctgs_li=f'{wkdir}/full_viral_ctgs.list'
    vs2_partial_ctgs_li=f'{wkdir}/vs2_partial_viral_ctgs.list'
    vb_partial_ctgs_li=f'{wkdir}/vb_partial_viral_ctgs.list'
    listToFile(full_ctgs,full_ctgs_li)
    listToFile(vs2_partial_ctgs,vs2_partial_ctgs_li)
    listToFile(vb_partial_ctgs,vb_partial_ctgs_li)
    return 0

if __name__=='__main__':
    if len(sys.argv)<2:
        print(f'{sys.argv[0]} <viral_identify_wkdir>')
    else:
        ctgList(sys.argv[1])
