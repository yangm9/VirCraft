#!/usr/bin/env python3

import sys,os
from ..config import setVari,conf
from ..process import cmdExec, general
from ..votus.deRep import VirRef 

class GeneFunc(VirRef):
    '''
    '''
    def __init__(self, config, outdir):
        VirRef.__init__(self, config, outdir)
        self.datadir = self.wkdir
        self.wkdir = f'{self.outdir}/06.func_annot'
        general.mkdir(self.wkdir)
        self.orfs_faa = f'{self.wkdir}/0.prodigal/all_votus.faa'
    def genePred(self):
        '''
        Predict genes from fasta file;
        2) Annotate the function of these genes.
        '''
        cmd = [self.envs]
        wkdir = f'{self.wkdir}/0.prodigal'
        general.mkdir(wkdir)
        temp_orfs_faa = f'{wkdir}/temp.orfs.faa'
        orfs_faa = f'{wkdir)/all_votus.faa'
        orfs_ffn = f'{wkdir}/all_votus.ffn'
        temp_orfs_ffn = f'{wkdir}/temp.orfs.ffn'
        temp_txt = f'{wkdir}/temp.txt'
        cmd.extend(
            ['prodigal', '-i', self.votus, '-a', temp_orfs_faa, '-d', temp_orfs_ffn, '-m', '-o', temp_txt, '-p meta -q\n', 
            'cut -f 1 -d \" \" ', temp_orfs_faa, '>', orfs_faa, '\n',
            'cut -f 1 -d \" \" ', temp_orfs_ffn, '>', orfs_ffn, '\n',
            f'rm -f {wkdir}/temp.*\n']
        )
        shell = f'{self.wkdir}/0.gene_predict.sh'
        general.printSH(shell, cmd)
        results = cmdExec.execute(cmd)
        return results
    def eggnogAnno(self):
        wkdir = f'{self.wkdir}/1.eggnog'
        general.mkdir(wkdir)
        cmd = [self.envs]
        anno_prefix = f'{wkdir}/all_votus'
        seed_orth = f'{anno_prefix}.emapper.seed_orthologs'
        eggout = f'{wkdir}/eggout'
        eggnog_db = self.confDict['EggNOGDB']
        cmd.extend(
            ['emapper.py -m diamond --no_annot --no_file_comments --cpu 40',
            '-i', orfs_faa, '-o', anno_prefix,
            '--data_dir', eggnog_db, '\n',
            'emapper.py', '--annotate_hits_table', seed_orth,
            '--no_file_comments', '-o', eggout, '--cpu 40', 
            '--data_dir', eggnog_db, '--override\n']        ]
        )
        shell = f'{wkdir}/eggnog_anno.sh'
        general.printSH(shell, cmd)
        results = cmdExec.execute(cmd)
        return results
    def keggAnno(self):
        wkdir = f'{func_anno_dir}/2.kegg'
        general.mkdir(wkdir)
        cmd = [self.envs]
        ko_prof = f'{self.confDict["kofamscanDB"]}/profiles'
        ko_list = f'{self.confDict["kofamscanDB"]}/ko_list'
        exec_anno = f'{wkdir}/all_votus.exec_annotation.txt'
        exec_anno_detail = f'{exec_anno}.xls'
        cmd.extend(
            ['exec_annotation -f mapper --cpu 16', '-p', ko_prof, '-k', ko_list,
            '-o', exec_anno, orfs_faa, '\n', 
            'exec_annotation -f detail --cpu 32', '-p', ko_prof, '-k', ko_list,
            '-o', exec_anno_detail, orfs_faa, '\n']
        )
        shell = f'{wkdir}/kegg_anno.sh'
        general.printSH(shell, cmd)
        results = cmdExec.execute(cmd)
        return results
