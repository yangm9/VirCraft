#!/usr/bin/env python3
#author: yangm@idsse.com
import os
import sys
import pandas as pd
import linkTab

FiltCondi="vs2_max_score>=0.9 or vb_prediction=='virus' or (dvf_score>=0.9 and dvf_pvalue<=0.1)"

PathDict={
    'virsorter2':'vs2-pass1/final-viral-score.tsv',
    'vibrant':'VIBRANT_{0}/VIBRANT_results_{0}/VIBRANT_machine_{0}.tsv',
    'deepvirfinder':'deepvirfinder'
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

CsvDict={
    'virsorter2':'{}/vs2_viral_cfgs.xls',    
    'vibrant':'{}/vb_viral_cfgs.xls',
    'deepvirfinder':'{}/dvf_viral_cfgs.xls'
}

ColsDict={
    'virsorter2':{
        'dsDNAphage':'vs2_dsDNAphage','ssDNA':'vs2_ssDNA','NCLDV':'vs2_NCLDV',
        'RNA':'vs2_RNA','lavidaviridae':'vs2_lavidaviridae',
        'max_score':'vs2_max_score','max_score_group':'vs2_max_score_group',
        'hallmark':'vs2_hallmark','viral':'vs2_viral',
        'cellular':'vs2_cellular'
    },
    'vibrant':{'prediction':'vb_prediction'},
    'deepvirfinder':{'score':'dvf_score','pvalue':'dvf_pvalue'}
}

#Get the fullpath results table from deepvirfinder or vibrant or virsorter2
def resultFile(name,tool,wkdir):
    wkdir=wkdir.rstrip('/')
    #name=os.path.basename(wkdir)
    result=wkdir+'/'+PathDict[tool]
    if tool=='deepvirfinder':
        file_name=os.listdir(result)[0]
        result+='/'+file_name
    elif tool=='vibrant':
        result=result.format(name)
    return result

def vCtgMerge(name,wkdir):
    full_ctgs=[]
    vs2_partial_ctgs=[]
    vb_partial_ctgs=[]
    for tool in PathDict.keys():
        result=resultFile(name,tool,wkdir)
        df=pd.read_csv(result,sep='\t')
        if tool=='virsorter2':
            df['seqname']=df['seqname'].apply(lambda x:x.rsplit('||',1)[0] if x.endswith('full') or x.endswith('lt2gene') else x)
            full_tmps=df.query('not seqname.str.endswith("partial")')['seqname'].tolist()
            partial_tmps=df.query('seqname.str.endswith("partial")')['seqname'].tolist()
            full_ctgs.extend(full_tmps)
            vs2_partial_ctgs.extend(partial_tmps)
            df.drop(columns=['length'], inplace=True)
        elif tool=='vibrant':
            full_tmps=df.query('not scaffold.str.contains("fragment")')['scaffold']
            partial_tmps=df.query('scaffold.str.contains("fragment")')['scaffold']
            full_ctgs.extend(full_tmps)
            vb_partial_ctgs.extend(partial_tmps)
        else:
            full_ctgs.extend(df[NameDict[tool]].tolist())
            df.drop(columns=['len'], inplace=True)
        df.rename(columns={NameDict[tool]:'Contig'},inplace=True)
        df.rename(columns=ColsDict[tool],inplace=True)
        csv_name=CsvDict[tool].format(wkdir)
        df.to_csv(csv_name,index=False,sep='\t')
    full_ctgs=list(set(full_ctgs))
    return full_ctgs,vs2_partial_ctgs,vb_partial_ctgs

def listToFile(list_l,list_f):
    LIST=open(list_f,'w')
    for ctg in list_l:
        LIST.write(f'{ctg}\n')
    return 0

def filtCtgList(all_merged_ctgs,filt_type):
    df=pd.read_csv(all_merged_ctgs,sep='\t')
    if filt_type=='tools':
        df['vs2_score']=df['vs2_max_score'].apply(lambda x:2 if x>=0.9 else (1 if x>0.7 else 0))
        df['vb_score']=df['vb_prediction'].apply(lambda x:1 if x=='virus' else 0)
        df['dvf_score']=df.apply(lambda x:1 if x['dvf_score']>=0.9 and x['dvf_pvalue']<=0.1 else 0, axis=1)
        df['evidences']=df['vs2_score']+df['vb_score']+df['dvf_score']
    elif(filt_type=='cutoff'):
        df=df.query(FiltCondi)
    else:
        pass
    postfix=f'.{filt_type}.xls'
    all_filt_ctgs=all_merged_ctgs.replace('.xls',postfix)
    df.to_csv(all_filt_ctgs,index=False,sep='\t')
    return 0

def ctgList(name,wkdir):
    full_ctgs,vs2_partial_ctgs,vb_partial_ctgs=vCtgMerge(name,wkdir)
    all_ctgs=['Contig']+full_ctgs+vs2_partial_ctgs+vb_partial_ctgs
    all_ctgs_li=f'{wkdir}/all_viral_ctgs.list'
    full_ctgs_li=f'{wkdir}/full_viral_ctgs.list'
    vs2_partial_ctgs_li=f'{wkdir}/vs2_partial_viral_ctgs.list'
    vb_partial_ctgs_li=f'{wkdir}/vb_partial_viral_ctgs.list'
    listToFile(all_ctgs,all_ctgs_li)
    listToFile(full_ctgs,full_ctgs_li)
    listToFile(vs2_partial_ctgs,vs2_partial_ctgs_li)
    listToFile(vb_partial_ctgs,vb_partial_ctgs_li)
    mark=all_ctgs_li
    for tool in CsvDict.keys():
        csv_name=CsvDict[tool].format(wkdir)
        merged_name=mark+'_'+tool
        linkTab.merge(mark,csv_name,'left','Contig',merged_name)
        os.remove(mark)
        mark=merged_name
    all_merged_ctgs=f'{wkdir}/all_viral_cfgs.xls'
    os.rename(mark,all_merged_ctgs)
    filtCtgList(all_merged_ctgs,'tools')
    filtCtgList(all_merged_ctgs,'cutoff')
    return 0

if __name__=='__main__':
    if len(sys.argv)<2:
        print(f'{sys.argv[0]} <viral_identify_wkdir>')
    else:
        ctgList(sys.argv[1],sys.argv[2])
