#!/usr/bin/env python3

import sys,os
from ..config import setVari,conf
from ..process import cmdExec, general
from ..fastqc import mergeRead

def funcAnno(config: str, outdir: str):
    '''
    1) Predict genes from fasta file;
    2) Annotate the function of these genes.
    '''
    groups, confDict, sampDict = conf.prepInfo(config)
    envs = setVari.selectENV('VirCraft')
    func_anno_cmd = [envs]
    func_anno_dir = f'{outdir}/06.func_annot'
    general.mkdir(func_anno_dir)
    fasta = f'{outdir}/03.vOTUs/merged_virus_positive_nodup.fa'
    prodigal_dir = f'{func_anno_dir}/0.prodigal'
    temp_orfs_faa = f'{prodigal_dir}/temp.orfs.faa'
    orfs_faa = f'{prodigal_dir}/merged_virus_positive_nodup.faa'
    orfs_ffn = f'{prodigal_dir}/merged_virus_positive_nodup.ffn'
    temp_orfs_ffn = f'{prodigal_dir}/temp.orfs.ffn'
    temp_txt = f'{prodigal_dir}/temp.txt'
    general.mkdir(prodigal_dir)
    func_anno_cmd.extend(
        ['prodigal', '-i', fasta, '-a', temp_orfs_faa, '-d', temp_orfs_ffn, '-m', '-o', temp_txt, '-p meta -q\n', 
        'cut -f 1 -d \" \" ', temp_orfs_faa, '>', orfs_faa, '\n',
        'cut -f 1 -d \" \" ', temp_orfs_ffn, '>', orfs_ffn, '\n',
        f'rm -f {prodigal_dir}/temp.*\n']
    )
    eggnog_dir = f'{func_anno_dir}/1.eggnog'
    general.mkdir(eggnog_dir)
    eggnog_anno_prefix = f'{eggnog_dir}/merged_virus_positive_nodup'
    seed_orth = f'{eggnog_dir}/merged_virus_positive_nodup.emapper.seed_orthologs'
    eggout = f'{eggnog_dir}/eggout'
    func_anno_cmd.extend(
        ['emapper.py -m diamond --no_annot --no_file_comments --cpu 40',
        '-i', orfs_faa, '-o', eggnog_anno_prefix,
        '--data_dir', confDict['EggNOGDB'], '\n',
        'emapper.py', '--annotate_hits_table', seed_orth,
        '--no_file_comments', '-o', eggout, '--cpu 40', 
        '--data_dir', confDict['EggNOGDB'], '--override\n']        ]
    )
    kegg_dir = f'{func_anno_dir}/2.kegg'
    general.mkdir(kegg_dir)
    ko_prof = f'{confDict['kofamscanDB']}/profiles'
    ko_list = f'{confDict['kofamscanDB']}/ko_list'
    exec_anno = f'{kegg_dir}/merged_virus_positive_nodup.exec_annotation.txt'
    exec_anno_detail = f'{exec_anno}.xls'
    func_anno_cmd.extend(
        ['exec_annotation -f mapper --cpu 16', '-p', ko_prof, '-k', ko_list,
        '-o', exec_anno, orfs_faa, '\n', 
        'exec_annotation -f detail --cpu 32', '-p', ko_prof, '-k', ko_list,
        '-o', exec_anno_detail, orfs_faa, '\n']
    )
    func_anno_sh = f'{func_anno_dir}/func_anno.sh'
    general.printSH(spades_sh, spades_cmd)
    results = cmdExec.execute(spades_cmd)
    return results
