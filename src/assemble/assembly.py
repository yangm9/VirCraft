#!/usr/bin/env python3

import os
import sys
from ..fastqc.readsQC import Reads 
from ..general import cmdExec,general

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    def __init__(self,fq1='',fq2='',outdir=''):
        super().__init__(fq1,fq2,outdir)
    def spades(self):
        '''
        Assemble metagenome by SPAdes.
        '''
        wkdir=f'{self.outdir}/spades'
        general.mkdir(wkdir)
        cmd=['spades.py','--pe1-1',self.fastqs[0],'--pe1-2',self.fastqs[1],
            '--careful','-t 32 -m 1300 -k 21,33,55,77,99,127','-o',wkdir,'\n']
        return cmd
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
            '-t','32','-m','80000000000','--tmp-dir',
            tmpdir,other_paras]
        return cmd
    def unmapReads(self):
        '''
        Align the FastQs back to Assembled Contigs.
        '''
        wkdir=f'{self.outdir}/alignment'
        general.mkdir(wkdir)
        scaffolds=f'{self.outdir}/spades/scaffolds.fasta'
        bwa_idx=f'{wkdir}/scaffoldsIDX'
        unused_sam=f'{wkdir}/unused_reads.sam'
        unused_fq=f'{wkdir}/unused_reads.fq'
        cmd=['bwa index -a bwtsw',scaffolds,'-p',bwa_idx,'\n',
            'bwa mem','-t 28',bwa_idx,self.fastqs[0],self.fastqs[1],
            '|grep -v NM:i:>',unused_sam,'\n',
            'sam_to_fastq.py',unused_sam,'>',unused_fq,'\n']
        return cmd,unused_fq
    def mergeFastA(self):
        scaffolds1=f'{self.outdir}/spades/scaffolds.fasta'
        scaffolds2=f'{self.outdir}/megahit/final.contigs.fa'
        scaffolds=f'{self.outdir}/final_assembly.fasta'
        stat_tab=f'{self.outdir}/stat.tab'
        cmd=['cat',scaffolds1,scaffolds2,'>',scaffolds,'\n',
             'assemb_stat.pl',scaffolds,scaffolds,f'>{stat_tab}\n']
        return cmd
    def filtFastA(self,cutoff=2000):
        '''
        Filter the fasta sequence by length (cutoff).
        '''
        wkdir=f'{self.outdir}/filter'
        scaffolds=f'{self.outdir}/final_assembly.fasta'
        filt_fa_prifix=f'{wkdir}/scaffolds.filt'
        cmd=['SeqLenCutoff.pl',scaffolds,filt_fa_prifix,str(cutoff),'\n']
        return cmd
    def Assemble(self):
        cmd=[self.envs]
        cmd.extend(self.spades())
        tmp_cmd,unused_fq=self.unmapReads()
        cmd.extend(tmp_cmd)
        cmd.extend(self.megahit([unused_fq]))
        cmd.extend(self.mergeFastA())
        cmd.extend(self.filtFastA(2000))
        cmd.extend(self.filtFastA(5000))
        cmd.extend(self.filtFastA(10000))
        shell=f'{self.outdir}/reads_assembly.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
