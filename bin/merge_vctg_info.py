#!/usr/bin/env python3
#author:   yangm@idsse.ac.cn
#versions: 0.0.1
#versions: 0.0.2 2023-07-26 22:42

# If you need to add new virus sequence identification software, the following content must be modified:
# 1. Add new key-value pairs to all dictionaries before the first function;
# 2. Modify the result file path in the resultFile function to align with the dictionaries;
# 3. Add the processing of results from the new software in the if...elif... structure within the vCtgMerge function;
# 4. Add a line in the calcCtgScore function to convert scores for the new software.

import os
import sys
import re
import pandas as pd
import linkTab

Tools = ['virsorter2', 'vibrant', 'genomad', 'deepvirfinder']

# The relative path to the final result table for viral identify tool, and it need to transform to the full path in this program.
PathDict = {
    'virsorter2' : 'vs2-pass1/final-viral-score.tsv',
    'vibrant' : 'VIBRANT_{0}/VIBRANT_results_{0}/VIBRANT_genome_quality_{0}.tsv',
    'genomad' : 'genomad/{0}_summary/{0}_virus_summary.tsv',
    'deepvirfinder' : 'deepvirfinder'
}

# The first columns name (Contig name) for final result table of each tool.
NameDict = {
    'virsorter2' : 'seqname',
    'vibrant' : 'scaffold',
    'deepvirfinder' : 'name',
    'genomad' : 'seq_name'
}

# The intermediate output of this program for each tool.
CsvDict = {
    'virsorter2' : '{}/vs2_viral_info.tsv',     
    'vibrant' : '{}/vb_viral_info.tsv',
    'deepvirfinder' : '{}/dvf_viral_info.tsv',
    'genomad' : '{}/gn_viral_info.tsv'
}

# The column dictionary for each tool,  and this will be used to rename the column names for the file named "all_viral_ctgs.score.tsv : ".
ColsDict = {
    'virsorter2' : {
        'dsDNAphage' : 'vs2_dsDNAphage', 'ssDNA' : 'vs2_ssDNA', 'NCLDV' : 'vs2_NCLDV',  'RNA' : 'vs2_RNA','lavidaviridae' : 'vs2_lavidaviridae', 
        'max_score' : 'vs2_max_score', 'max_score_group' : 'vs2_max_score_group', 'hallmark' : 'vs2_hallmark', 'viral' : 'vs2_viral', 'cellular' : 'vs2_cellular'
    },
    'vibrant' :  {'type' : 'vb_lifestyle'},
    'deepvirfinder' : {'score' : 'dvf_v_score', 'pvalue' : 'dvf_pvalue'},
    'genomad' :  {
        'topology' : 'gn_topology', 'coordinates' : 'gn_coordinates', 'n_genes' : 'gn_genes', 'genetic_code' : 'gn_genetic_code', 'virus_score' : 'gn_v_score',
        'fdr' : 'gn_fdr', 'n_hallmarks' : 'gn_hallmarks', 'marker_enrichment' : 'gn_marker_enrichment', 'taxonomy' : 'gn_taxonomy'
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
    result = wkdir + '/' + PathDict[tool]
    if tool == 'deepvirfinder':
        file_name = os.listdir(result)[0]
        result += '/' + file_name
    elif tool == 'vibrant' or tool == 'genomad':
        result = result.format(name)
    return result

#Extract a list of full and partial contigs from the results of dvf, vb, vs2, and genomad tools, meanwhile output key information for each tools
def vCtgMerge(name, wkdir):
    all_ctgs = []
    for tool in Tools:
        result = resultFile(name, tool, wkdir)
        df = pd.read_csv(result, sep='\t', header=0)
        if tool == 'virsorter2':
            df[['seqname', 'vs2_partial']] = df['seqname'].str.split(r'\|\|', expand=True)
            df.drop(columns=['length'], inplace=True)
        elif tool == 'vibrant':
            df = df[['scaffold', 'type']].drop_duplicates()
            df['vb_isPhage'] = 1
            try: #In case the VIBRANT doesn't output fragment results
                df[['scaffold', 'vb_partial']] = df['scaffold'].str.split(r'_frag', expand=True)
                df['vb_partial'] = df['vb_partial'].str.replace('ment_', 'fragment_')
            except ValueError:
                df['vb_partial'] = ''

            agg_dict = {
                'type': 'first',
                'vb_isPhage': 'max',  # Take the maximum value of vb_isPhage
                'vb_partial': lambda x: ','.join(str(item) for item in x if item is not None)
            }
            df = df.groupby('scaffold', as_index=False).agg(agg_dict)
        elif tool == 'deepvirfinder':
            df.drop(columns = ['len'], inplace=True)
        elif tool == 'genomad':
            try:
                df[['seq_name', 'gn_partial']] = df['seq_name'].str.split(r'\|pro', expand=True)
                df['gn_partial'] = df['gn_partial'].str.replace('virus_', 'provirus_')
            except ValueError:
                df['gn_partial'] = ''
            df.drop(columns = ['length'], inplace=True)
        df.rename(columns = {NameDict[tool]: 'Contig'}, inplace=True)
        all_ctgs.extend(df['Contig'].tolist())
        all_ctgs = list(set(all_ctgs))
        df.rename(columns=ColsDict[tool], inplace=True)
        csv_name = CsvDict[tool].format(wkdir)
        df.to_csv(csv_name, index=False, sep='\t')
    return all_ctgs

#Calculate the final score according to all results.
def calcCtgScore(all_merged_ctg_tsv, filt_mode='permissive'):
    df = pd.read_csv(all_merged_ctg_tsv, sep='\t')
    df['vb_lifestyle'].fillna('Undetermined', inplace=True)
    df['gn_taxonomy'].fillna('Unclassified', inplace=True)
    df['vs2_score'] = df.apply(lambda x : 2 if x['vs2_max_score'] >= 0.9 else (1 if x['vs2_max_score'] >= 0.7 and x['vs2_hallmark'] >= 1 else 0), axis = 1)
    df['vb_score'] = df['vb_isPhage'].apply(lambda x : 1 if x == 1 else 0)
    df['dvf_score'] = df.apply(lambda x : 1 if x['dvf_v_score'] >= 0.9 and x['dvf_pvalue'] <= 0.05 else 0, axis = 1)
    df['gn_score'] = df.apply(lambda x : 2 if x['gn_v_score'] >= 0.8 else (1 if x['gn_v_score'] >= 0.7 and x['gn_hallmarks'] >= 1 else 0), axis = 1)
    df['score'] = df['vs2_score'] + df['vb_score'] + df['dvf_score'] + df['gn_score']
    postfix = '.score.tsv'
    all_ctg_score_tsv = all_merged_ctg_tsv.replace('.tsv', postfix)
    df.to_csv(all_ctg_score_tsv, index=False, sep='\t')
    score_cutoff_dict = {'permissive': 1, 'strict': 2}
    postfix = f'.score.{filt_mode}.tsv'
    all_ctg_score_filt_tsv = all_merged_ctg_tsv.replace('.tsv', postfix)
    df = df.query(f'score >= {score_cutoff_dict[filt_mode]}')
    df.to_csv(all_ctg_score_filt_tsv, index=False, sep='\t')
    return 0

#Main function
def ctgList(name, wkdir, filt_mode='permissive'):
    all_nh_ctgs = vCtgMerge(name, wkdir)
    all_ctgs = ['Contig'] + all_nh_ctgs # the 1st column of output
    all_ctgs_li = f'{wkdir}/all_viral_ctgs.list'
    listToFile(all_ctgs, all_ctgs_li)
    mark = all_ctgs_li
    for tool in Tools:
        csv_name = CsvDict[tool].format(wkdir) #{vb}_viral_cfgs.tsv
        merged_name = mark + '_' + tool # all_viral_ctgs.list_virsorter2_vibrant_deepvirfinder_genomad 
        linkTab.merge(mark, csv_name, 'left', 'Contig', merged_name)
        if mark != all_ctgs_li: # To prevent accidental deletion of all_viral_ctgs.list
            os.remove(mark)
        mark = merged_name
    all_merged_ctg_tsv = f'{wkdir}/all_viral_ctgs.tsv'
    os.rename(mark, all_merged_ctg_tsv)
    calcCtgScore(all_merged_ctg_tsv, filt_mode)
    #os.remove(all_merged_ctg_tsv)
    return 0

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'{sys.argv[0]} <metagenomic_name> <viral_identify_wkdir> [filter_mode]')
    elif len(sys.argv) == 2:
        ctgList(sys.argv[1], sys.argv[2])
    else:
        ctgList(sys.argv[1], sys.argv[2], sys.argv[3])
