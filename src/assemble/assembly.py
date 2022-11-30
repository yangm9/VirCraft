#!/usr/bin/env python3

import os
import sys
from ..general import cmdExec,general
from ..fastqc.reads import Reads

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,*args,config,outdir,**kwargs):
        Reads.__init__(self,config,outdir)
        self.wkdir=f'{self.outdir}/01.assembly'
    def spades(self,fastqs:list,group:str):
        '''
        Assemble metagenome by SPAdes for single group.
        '''
        cmd=[self.envs]
        wkdir=f'{self.wkdir}/{group}'
        general.mkdir(wkdir)
        cmd.extend(
            ['spades.py','--pe1-1',fastqs[0],'--pe1-2',fastqs[1],
             '--careful','-t 30 -m 1300 -k 21,33,55,77,99,127',
             '-o',wkdir]
        )
        shell=f'{self.wkdir}/{group}_spades.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def megahit(self,fastqs:list,group:str):
        '''
        Assemble metagenome by SPAdes for single group.
        '''
        cmd=[self.envs]
        input_para=''
        other_paras=''
        wkdir=f'{self.wkdir}/{group}'
        tmpdir=f'{wkdir}/megahit.tmp'
        outdir=f'{wkdir}/megahit'
        if len(fastqs)==1:
            input_para=f'-r {fastqs[0]}'
        else:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]}'
            other_paras='--continue'
        cmd.extend(
            ['megahit',input_para,'-o',outdir,
             '-t','32','-m','80000000000','--tmp-dir',tmpdir,other_paras]
        )
        shell=f'{self.wkdir}/{group}_megahit.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def filtFastA(self,grp:str,cutoff:int):
        '''
        Filter the fasta sequence by length (cutoff).
        '''
        wkdir=f'{self.wkdir}/{grp}'
        scaffolds=f'{wkdir}/scaffolds.fasta'
        filt_fa_prifix=f'{wkdir}/scaffolds.filt'
        filt_cmd=['SeqLenCutoff.pl',scaffolds,filt_fa_prifix,cutoff]
        filt_sh=f'{self.wkdir}/{grp}_filt_scaffolds.sh'
        general.printSH(filt_sh,filt_cmd)
        results=cmdExec.execute(filt_cmd)
        return results
    def statFastA(self,grp:str):
        wkdir=f'{self.wkdir}/{grp}'
        contigs=f'{wkdir}/scaffolds.fasta'
        scaffolds=f'{wkdir}/scaffolds.fasta'
        stat_tab=f'{wkdir}/stat.tab'
        cmd=['assemb_stat.pl',contigs,scaffolds,f'>{stat_tab}\n']
        shell=f'{self.wkdir}/{grp}_fasta_stat.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def unmapReads(self,fastqs:list,grp:str):
        '''
        Align the FastQs back to Assembled Contigs.
        '''
        cmd=[self.envs]
        wkdir=f'{self.wkdir}/{grp}'
        scaffolds=f'{wkdir}/scaffolds.fasta'
        bwa_idx=f'{wkdir}/scaffoldsIDX'
        unused_sam=f'{wkdir}/unused_by_spades_{grp}.sam'
        unused_fq=f'{wkdir}/unused_by_spades_{grp}.fq'
        cmd.extend(
            ['bwa index -a bwtsw',scaffolds,'-p',bwa_idx,'\n',
             'bwa mem','-t 28',bwa_idx,fastqs[0],fastqs[1],
             '|grep -v NM:i:>',unused_sam,'\n',
             'sam_to_fastq.py',unused_sam,'>',unused_fq,'\n']
        )
        shell=f'{self.wkdir}/{grp}_unmapped_reads.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return unused_fq
    def Assemble(self):
        results=''
        for grp in self.groups:
            fastq_1=f'{self.fq_dir}/{grp}_1.fq'
            fastq_2=f'{self.fq_dir}/{grp}_2.fq'
            fastqs=[fastq_1,fastq_2]
            results+=self.spades(fastqs,grp)
            results+=self.unmapReads(fastqs,grp)
            results+=self.megahit(fastqs,grp)
            results+=self.filtFastA(grp,'2000')
            results+=self.filtFastA(grp,'5000')
            results+=self.filtFastA(grp,'10000')
            results+=self.statFastA(grp)
        return results
