#!/usr/bin/env python3
#author:   yangm@idsse.com
#versions: 0.0.1
#versions: 0.0.2 2023-07-26 22:42

import os
import sys
import re
import pandas as pd
import linkTab

#The relative path to the final result table for viral identify tool, and it need to transform to the full path in this program.
PathDict={
    'virsorter2':'vs2-pass1/final-viral-score.tsv',
    'vibrant':'VIBRANT_{0}/VIBRANT_results_{0}/VIBRANT_machine_{0}.tsv',
    'deepvirfinder':'deepvirfinder'
}

#The first columns name (Contig name) for final result table of each tool.
NameDict={
    'virsorter2':'seqname',
    'vibrant':'scaffold',
    'deepvirfinder':'name'
}

#The Final output FastA file from each tool.
FastaDict={
    'virsorter2':'vs2-pass1/final-viral-combined.fa',
    'vibrant':'VIBRANT_{0}/VIBRANT_phages_{0}/{0}.phages_combined.fna',
    'deepvirfinder':''
}

#The intermediate output from this program for each tool.
CsvDict={
    'virsorter2':'{}/vs2_viral_info.tsv',    
    'vibrant':'{}/vb_viral_info.tsv',
    'deepvirfinder':'{}/dvf_viral_info.tsv'
}

#The column dictionary for each tool, and this will be used to rename the column names for the file named "all_viral_ctgs.score.tsv:".
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

#Output a list to a file
def listToFile(list_l,list_f):
    LIST=open(list_f,'w')
    for ctg in list_l:
        LIST.write(f'{ctg}\n')
    return 0

#Get the fullpath results table from deepvirfinder, vibrant or virsorter2
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

#Extract a list of full and partial contigs from the results of dvf, vb and vs2 tools, meanwhile output key information for each tools.
def vCtgMerge(name,wkdir):
    all_ctgs=[]
    for tool in PathDict.keys():
        result=resultFile(name,tool,wkdir)
        df=pd.read_csv(result,sep='\t')
        if tool=='virsorter2':
            df[['seqname','vs2_partial']]=df['seqname'].str.split(r'\|\|',expand=True)
            df.drop(columns=['length'], inplace=True)
        elif tool=='vibrant':
            phage_txt=FastaDict[tool].format(name).replace('.fna','.txt')
            phage_txt=f'{wkdir}/{phage_txt}'
            phage_df=pd.read_csv(phage_txt,sep='\t',header=None)
            phage_list=phage_df[0].tolist()
            df['vb_isPhage']=df['scaffold'].isin(phage_list).astype(int)
            try: #In case the VIBRANT doesn't output fragment results
                df[['scaffold','vb_partial']]=df['scaffold'].str.split(r'_frag',expand=True)
                df['vb_partial']=df['vb_partial'].str.replace('ment_','fragment_')
            except ValueError:
                df['vb_partial'] = ''
        else:
            df.drop(columns=['len'], inplace=True)
        df.rename(columns={NameDict[tool]:'Contig'},inplace=True)
        all_ctgs.extend(df['Contig'].tolist())
        all_ctgs=list(set(all_ctgs))
        df.rename(columns=ColsDict[tool],inplace=True)
        csv_name=CsvDict[tool].format(wkdir)
        df.to_csv(csv_name,index=False,sep='\t')
    return all_ctgs

#Calculate the final score according to all results.
def calcCtgScore(all_merged_ctgs):
    df=pd.read_csv(all_merged_ctgs,sep='\t')
    df['vs2_score']=df['vs2_max_score'].apply(lambda x:2 if x>=0.9 else (1 if x>=0.7 else 0))
    df['vb_score']=df['vb_isPhage'].apply(lambda x:1 if x==1 else 0)
    df['dvf_scores']=df.apply(lambda x:1 if x['dvf_score']>=0.9 and x['dvf_pvalue']<=0.1 else 0, axis=1)
    df['score']=df['vs2_score']+df['vb_score']+df['dvf_scores']
    postfix=f'.score.tsv'
    all_filt_ctgs=all_merged_ctgs.replace('.tsv',postfix)
    df.to_csv(all_filt_ctgs,index=False,sep='\t')
    return 0

#Main Function
def ctgList(name,wkdir):
    all_nh_ctgs=vCtgMerge(name,wkdir)
    all_ctgs=['Contig']+all_nh_ctgs
    all_ctgs_li=f'{wkdir}/all_viral_ctgs.list'
    listToFile(all_ctgs,all_ctgs_li)
    mark=all_ctgs_li
    for tool in CsvDict.keys():
        csv_name=CsvDict[tool].format(wkdir) #{vb}_viral_cfgs.tsv
        merged_name=mark+'_'+tool
        linkTab.merge(mark,csv_name,'left','Contig',merged_name)
        if mark!=all_ctgs_li: os.remove(mark)
        mark=merged_name
    all_merged_ctgs=f'{wkdir}/all_viral_ctgs.tsv'
    os.rename(mark,all_merged_ctgs)
    calcCtgScore(all_merged_ctgs)
    return 0

if __name__=='__main__':
    if len(sys.argv)<2:
        print(f'{sys.argv[0]} <metagenomic_name> <viral_identify_wkdir>')
    else:
        ctgList(sys.argv[1],sys.argv[2])
