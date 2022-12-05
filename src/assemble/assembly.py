#!/usr/bin/env python3

import os
import sys
from ..general import cmdExec,general
from ..fastqc.readsQC import Reads 

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(threads)
        self.methDict={
            's':self.spades(),
            'm':self.megahit()
        }
    def spades(self,fastqs:list):
        '''
        Assemble metagenome by SPAdes.
        '''
        wkdir=f'{self.outdir}/spades'
        general.mkdir(wkdir)
        cmd=['spades.py','--pe1-1',fastqs[0],'--pe1-2',fastqs[1],
            '-t',self.threads,'-o',wkdir,
            '--careful -m 1300 -k 21,33,55,77,99,127','\n']
        scaf=f'{wkdir}/scaffolds.fasta'
        return cmd,scaf
    def megahit(self,fastqs:list):
        '''
        Assemble metagenome by megahit.
        '''
        input_para=''
        other_paras=''
        wkdir=f'{self.outdir}/megahit'
        general.mkdir(wkdir)
        tmpdir=f'{wkdir}/megahit.tmp'
        general.mkdir(tmpdir)
        if len(fastqs)==1:
            input_para=f'-r {fastqs[0]}'
        else:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]}'
            other_paras='--continue'
        cmd=['megahit',input_para,'-o',self.outdir,
            '-t',self.threads,'-m','80000000000','--tmp-dir',
            tmpdir,other_paras]
        scaf=f'{wkdir}/final.contigs.fa'
        return cmd,scaf
    def unmapReads(self,scaf:str):
        '''
        Align the FastQs back to Assembled Contigs.
        '''
        wkdir=f'{self.outdir}/alignment'
        general.mkdir(wkdir)
        bwa_idx=f'{wkdir}/scaffoldsIDX'
        unused_sam=f'{wkdir}/unused_reads.sam'
        unused_fq=f'{wkdir}/unused_reads.fq'
        cmd=['bwa index -a bwtsw',scaf,'-p',bwa_idx,'\n',
            'bwa mem','-t',self.threads,bwa_idx,self.fastqs[0],self.fastqs[1],
            '|grep -v NM:i:>',unused_sam,'\n',
            'sam_to_fastq.py',unused_sam,'>',unused_fq,'\n']
        return cmd,unused_fq
    def mergeFastA(self,scafs:list):
        scaf=f'{self.outdir}/final_assembly.fasta'
        stat_tab=f'{self.outdir}/stat.tab'
        cmd=['cat',scafs[0],scafs[1],'>',scaf,'\n',
             'assemb_stat.pl',scaffolds,scaffolds,f'>{stat_tab}\n']
        return cmd,scaf
    def filtFastA(self,cutoff=5000):
        '''
        Filter the fasta sequence by length (cutoff).
        '''
        wkdir=f'{self.outdir}/filter'
        scaf=f'{self.outdir}/final_assembly.fasta'
        filt_fa_prifix=f'{wkdir}/scaffolds.filt'
        cmd=['SeqLenCutoff.pl',scaf,filt_fa_prifix,str(cutoff),'\n']
        return cmd
    def mixAsse(self,fastqs,process='sm'):
        cmd=[]
        scafs=[]
        step_num=len(process)
        tmp_cmd,tmp_scaf=self.methDict[process[0]](fastqs)
        cmd.extend(tmp_cmd)
        scafs.append(tmp_scaf)
        if step_num==2:
            tmp_cmd,unused_fq=self.unmapReads(tmp_scaf)
            cmd.extend(tmp_cmd)
            tmp_cmd,tmp_scaf=self.methDict[process[1]]([unused_fq])
            cmd.extend(tmp_cmd)
        scafs.append(tmp_scaf)
        cmd.extend(self.mergeFastA(scafs))
        return cmd,scafs
    def Assemble(self,process='sm',cutoff=5000,threads=8):
        cmd=[self.envs]
        scafs=[]
        tmp_cmd,scafs=mixAsse(self.fastqs,process)
        cmd.extend(self.filtFastA(cutoff))
        shell=f'{self.outdir}/{self.samp}_assembly.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
