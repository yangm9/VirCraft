#!/usr/bin/env python3

import os
import sys
import glob
from ..config.config import Reads,Seq
from ..general import cmdExec,general
from ..identify.viridsop import VirScan

class AlnData(Reads):
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(threads)
    def alnReads(self,samp:str,bwa_idx:str):
        '''
        Align the reads to the vOTUs for each sample.
        '''
        cmd=[self.envs]
        raw_bam=f'{self.outdir}/{samp}.raw.bam'
        sort_bam=f'{self.outdir}/{samp}.sort.bam'
        cmd.extend(
            ['bwa mem','-t',self.threads,bwa_idx,
            self.fastqs[0],self.fastqs[1],
            '|samtools view','-o',raw_bam,'-@ 28 -b -S\n',
            'samtools sort',raw_bam,'-o',sort_bam,'-@ 28\n']
        )
        shell=f'{self.outdir}/{samp}_aln.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results

class multiAln(Seq):
    def __init__(self,fastqs='',outdir='',threads=8)
        super().__init__(fasta,outdir)
        self.fastqs=glob.glob(os.path.abspath(fastqs))
        self.indir=os.path.dirname(self.fastqs[0])
        self.postfix=os.path.splitext(self.fastqs[0])[1]
        fq1postfix=f'_1{self.postfix}'
        self.samps=[os.path.basename(name.replace(fq1postfix,'')) for name in fastqs if name.endswith(fq1postfix)]
        self.threads=str(threads)
    def AlnBySamp(self):
        results=''
        _,bwa_idx=self.mkBwaIdx
        for samp in self.samps:
            fq1=f'{self.indir}/{samp}_1{self.postfix}'
            fq2=f'{self.indir}/{samp}_2{self.postfix}'
            Aln=AlnData(fq1,fq2,self.outdir,self.threads)
            results+=Aln.alnReads(samp,bwa_idx)
        return results
