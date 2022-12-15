#!/usr/bin/env python3

import os
import sys
from ..general import cmdExec,general
from ..dataqc.fastqc import Reads 

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(threads)
        self.methDict={
            's':self.spades,
            'm':self.megahit
        }
    def spades(self,fastqs:list):
        '''
        Assemble metagenome by SPAdes.
        '''
        wkdir=f'{self.outdir}/spades'
        general.mkdir(wkdir)
        cmd=[]
        if len(fastqs)==1:
            cmd=['spades.py','-s',fastqs[0]]
        elif len(fastqs)==2:
            cmd=['spades.py','--pe1-1',fastqs[0],'--pe1-2',fastqs[1]]
        else:
            pass
        cmd.extend(
            ['-t',self.threads,'-o',wkdir,
            '--careful -m 1300 -k 21,33,55,77,99,127','\n']
        )
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
        elif len(fastqs)==2:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]}'
            other_paras='--continue'
        else:
            pass
        cmd=['megahit',input_para,'-o',self.outdir,
            '-t',self.threads,'-m','80000000000',
            '--tmp-dir',tmpdir,other_paras,'\n']
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
             'assemb_stat.pl',scaf,scaf,f'>{stat_tab}\n']
        return cmd,scaf
    def statFilt(self,scaf,cutoff=5000):
        '''
        Filter the fasta sequence by length (cutoff).
        '''
        wkdir=f'{self.outdir}/filter'
        general.mkdir(wkdir)
        stat_prefix=f'{wkdir}/scaffolds.stat'
        filt_prefix=f'{wkdir}/scaffolds.filt'
        cmd=['fasta_size_distribution_plot.py',scaf,'-o',stat_prefix,
            '-t Sequence Size Distribution -s 2000 -g 10\n',
            'SeqLenCutoff.pl',scaf,filt_prefix,str(cutoff),'\n'
        ]
        return cmd
    def mixAsse(self,fastqs,process='sm'):
        '''
        Analysis the Assembly process accordding to the process set.
        '''
        cmd=[]
        scafs=[]
        tmp_cmd,tmp_scaf=self.methDict[process[0]](fastqs)
        cmd.extend(tmp_cmd)
        scafs.append(tmp_scaf)
        if len(process)==2:
            tmp_cmd,unused_fq=self.unmapReads(tmp_scaf)
            cmd.extend(tmp_cmd)
            tmp_cmd,tmp_scaf=self.methDict[process[1]]([unused_fq])
            cmd.extend(tmp_cmd)
            scafs.append(tmp_scaf)
        tmp_cmd,scaf=self.mergeFastA(scafs)
        cmd.extend(tmp_cmd)
        return cmd,scaf
    def Assemble(self,process='sm',cutoff=5000):
        cmd=[self.envs]
        scafs=[]
        tmp_cmd,scaf=self.mixAsse(self.fastqs,process)
        cmd.extend(tmp_cmd)
        cmd.extend(self.statFilt(scaf,cutoff))
        shell=f'{self.outdir}/{self.samp}_assembly.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
