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

'''
如果要添加新的病毒序列鉴定软件，需要修改以下内容：
1. 第一个函数前所有的字典中需要添加新的键值对；
2. resultFile函数中的结果文件位置需要配合字典进行修改
3. vCtgMerge函数中if……elif……结构中需要对应添加新软件的结果处理
4. calcCtgScore函数中添加一行对新软件进行分数转换
'''

PathDict = {
    'virsorter2' : 'vs2-pass1/final-viral-score.tsv',
    'vibrant' : 'VIBRANT_{0}/VIBRANT_results_{0}/VIBRANT_machine_{0}.tsv',
    'genomad' : 'genomad/{0}_summary/{0}_virus_summary.tsv',
    'deepvirfinder' : 'deepvirfinder'
}

#The first columns name (Contig name) for final result table of each tool.
NameDict = {
    'virsorter2' : 'seqname',
    'vibrant' : 'scaffold',
    'deepvirfinder' : 'name',
    'genomad' : 'seq_name'
}

#The Final output FastA file from each tool.
FastaDict = {
    'virsorter2' : 'vs2-pass1/final-viral-combined.fa',
    'vibrant' : 'VIBRANT_{0}/VIBRANT_phages_{0}/{0}.phages_combined.fna',
    'deepvirfinder' : '',
    'genomad' : ''
}

#The intermediate output from this program for each tool.
CsvDict = {
    'virsorter2' : '{}/vs2_viral_info.tsv',     
    'vibrant' : '{}/vb_viral_info.tsv',
    'deepvirfinder' : '{}/dvf_viral_info.tsv',
    'genomad' : '{}/gm_viral_info.tsv'
}

#The column dictionary for each tool,  and this will be used to rename the column names for the file named "all_viral_ctgs.score.tsv : ".
ColsDict = {
    'virsorter2' : {
        'dsDNAphage' : 'vs2_dsDNAphage', 'ssDNA' : 'vs2_ssDNA', 'NCLDV' : 'vs2_NCLDV',  'RNA' : 'vs2_RNA','lavidaviridae' : 'vs2_lavidaviridae', 
        'max_score' : 'vs2_max_score', 'max_score_group' : 'vs2_max_score_group', 'hallmark' : 'vs2_hallmark', 'viral' : 'vs2_viral', 'cellular' : 'vs2_cellular'
    },
    'vibrant' :  {'prediction' : 'vb_prediction'},
    'deepvirfinder' : {'score' : 'dvf_v_score', 'pvalue' : 'dvf_pvalue'},
    'genomad' :  {
        'topology' : 'gm_topology', 'coordinates' : 'gm_coordinates', 'n_genes' : 'gm_genes', 'genetic_code' : 'gm_genetic_code', 'virus_score' : 'gm_v_score',
        'fdr' : 'gm_fdr', 'n_hallmarks' : 'gm_hallmarks', 'marker_enrichment' : 'gm_marker_enrichment', 'taxonomy' : 'gm_taxonomy'
    }
}

#Output a list to a file
def listToFile(list_l, list_f):
    LIST = open(list_f, 'w')
    for ctg in list_l:
        LIST.write(f'{ctg}\n')
    return 0

#Get the fullpath results table from deepvirfinder, vibrant, virsorter2 or genomad
def resultFile(name, tool, wkdir):
    wkdir = wkdir.rstrip('/')
    result = wkdir + '/'+PathDict[tool]
    if tool == 'deepvirfinder':
        file_name = os.listdir(result)[0]
        result += '/' + file_name
    elif tool == 'vibrant' or tool == 'genomad':
        result = result.format(name)
    return result

#Extract a list of full and partial contigs from the results of dvf, vb and vs2 tools, meanwhile output key information for each tools.
def vCtgMerge(name, wkdir):
    all_ctgs = []
    for tool in PathDict.keys():
        result = resultFile(name, tool, wkdir)
        df = pd.read_csv(result, sep='\t')
        if tool == 'virsorter2':
            df[['seqname', 'vs2_partial']] = df['seqname'].str.split(r'\|\|', expand=True)
            df.drop(columns=['length'], inplace=True)
        elif tool == 'vibrant':
            phage_txt = FastaDict[tool].format(name).replace('.fna', '.txt')
            phage_txt = f'{wkdir}/{phage_txt}'
            try:
                phage_df = pd.read_csv(phage_txt, sep='\t', header=None)
                phage_list = phage_df[0].tolist()
            except pd.errors.EmptyDataError:
                phage_df = pd.DataFrame()
                phage_list = []
            df['vb_isPhage'] = df['scaffold'].isin(phage_list).astype(int)
            try: #In case the VIBRANT doesn't output fragment results
                df[['scaffold', 'vb_partial']] = df['scaffold'].str.split(r'_frag', expand=True)
                df['vb_partial'] = df['vb_partial'].str.replace('ment_', 'fragment_')
            except ValueError:
                df['vb_partial'] = ''
        elif tool == 'deepvirfinder':
            df.drop(columns = ['len'], inplace = True)
        elif tool == 'genomad':
            try:
                df[['seq_name', 'gm_partial']] = df['seq_name'].str.split(r'\|pro', expand = True)
                df['gm_partial'] = df['gm_partial'].str.replace('virus_', 'provirus_')
            except ValueError:
                df['gm_partial'] = ''
            df.drop(columns = ['length'], inplace = True)
        df.rename(columns = {NameDict[tool]: 'Contig'}, inplace = True)
        all_ctgs.extend(df['Contig'].tolist())
        all_ctgs = list(set(all_ctgs))
        df.rename(columns = ColsDict[tool], inplace = True)
        csv_name = CsvDict[tool].format(wkdir)
        df.to_csv(csv_name, index = False, sep = '\t')
    return all_ctgs

#Calculate the final score according to all results.
def calcCtgScore(all_merged_ctgs):
    df=pd.read_csv(all_merged_ctgs, sep='\t')
    df['vs2_score'] = df.apply(lambda x : 2 if x['vs2_max_score'] >= 0.9 else (1 if x['vs2_max_score'] >= 0.7 and x['vs2_hallmark'] > 0 else 0), axis = 1)
    df['vb_score'] = df['vb_isPhage'].apply(lambda x : 1 if x == 1 else 0)
    df['dvf_score'] = df.apply(lambda x : 1 if x['dvf_v_score'] >= 0.9 and x['dvf_pvalue'] <= 0.1 else 0, axis = 1)
    df['gm_score'] = df.apply(lambda x : 2 if x['gm_v_score'] >= 0.8 else (1 if x['gm_v_score'] >= 0.7 and x['gm_hallmarks'] > 0 else 0), axis = 1)
    df['score'] = df['vs2_score'] + df['vb_score'] + df['dvf_score'] + df['gm_score']
    postfix = f'.score.tsv'
    all_filt_ctgs = all_merged_ctgs.replace('.tsv', postfix)
    df.to_csv(all_filt_ctgs, index=False, sep='\t')
    return 0

#Main function
def ctgList(name, wkdir):
    all_nh_ctgs = vCtgMerge(name, wkdir)
    all_ctgs = ['Contig'] + all_nh_ctgs
    all_ctgs_li = f'{wkdir}/all_viral_ctgs.list'
    listToFile(all_ctgs, all_ctgs_li)
    mark = all_ctgs_li
    for tool in CsvDict.keys():
        csv_name = CsvDict[tool].format(wkdir) #{vb}_viral_cfgs.tsv
        merged_name = mark + '_' + tool
        linkTab.merge(mark, csv_name, 'left', 'Contig', merged_name)
        if mark != all_ctgs_li:
            os.remove(mark)
        mark = merged_name
    all_merged_ctgs = f'{wkdir}/all_viral_ctgs.tsv'
    os.rename(mark, all_merged_ctgs)
    calcCtgScore(all_merged_ctgs)
    return 0

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'{sys.argv[0]} <metagenomic_name> <viral_identify_wkdir>')
    else:
        ctgList(sys.argv[1], sys.argv[2])
