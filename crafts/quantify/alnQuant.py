#!/usr/bin/env python3

import os
import sys
import glob
from ..config.config import Reads,Seq
from ..general import cmdExec,general
from ..identify.viridsop import VirScan

class VirCount(Reads):
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(threads)
    def aln(self,samp:str,bwa_idx:str):
        '''
        Align the reads to the vOTUs for each sample.
        '''
        raw_bam=f'{self.outdir}/{samp}.raw.bam'
        sort_bam=f'{self.outdir}/{samp}.sort.bam'
        cmd=['bwa mem','-t',self.threads,bwa_idx,
            self.fastqs[0],self.fastqs[1],
            '|samtools view','-o',raw_bam,'-@ 28 -b -S\n',
            'samtools sort',raw_bam,'-o',sort_bam,'-@ 28\n']
        return cmd
    def coverm(self,samp:str):
        tpm=f'{self.outdir}/{samp}.tpm'
        sort_bam=f'{self.outdir}/{samp}.sort.bam'
        cmd=['coverm contig','-b',sort_bam,'-t',self.threads,
            '--min-read-aligned-length 50 --min-read-percent-identity 0.95 --proper-pairs-only -m tpm',
            '>',tpm]
        return cmd