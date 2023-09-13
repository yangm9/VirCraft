#!/usr/bin/env python3

import os
import sys
import glob
from ..data.bioseq import Reads
from ..data.bioseq import Seq
from ..identify.viridsop import VirScan

class VirCount(Reads):
    coverm_args='--min-read-percent-identity 0.95 --min-read-aligned-length 50 --min-covered-fraction 10 -m mean'
    #coverm_args='--min-read-percent-identity 0.95 --min-read-aligned-length 50 --min-covered-fraction 10 --proper-pairs-only -m mean'
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(int(threads)//self.BATCH_SIZE)
    def bwa(self,samp:str,bwa_idx:str):
        '''
        Align the reads to the vOTUs for each sample.
        '''
        raw_bam=f'{self.outdir}/{samp}.raw.bam'
        sort_bam=f'{self.outdir}/{samp}.sort.bam'
        cmd=[utils.selectENV('VC-Quantify')]
        cmd.extend(
            ['bwa mem','-t',self.threads,bwa_idx,
            self.fastqs[0],self.fastqs[1],
            '|samtools view','-o',raw_bam,'-@ 28 -b -S\n',
            'samtools sort',raw_bam,'-o',sort_bam,'-@ 28\n',
            'samtools index',sort_bam,'\n']
        )
        return cmd
    def coverm(self,samp:str):
        cov=f'{self.outdir}/{samp}.cov'
        sort_bam=f'{self.outdir}/{samp}.sort.bam'
        cmd=[utils.selectENV('VC-Quantify')]
        cmd.extend(
            ['coverm contig','-b',sort_bam,'-t',self.threads,
            self.coverm_args,'>',cov,'\n']
        )
        return cmd

class GeneCount(Reads):
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(int(threads)//self.BATCH_SIZE)
    def salmon(self,samp:str,salmon_idx:str):
        wkdir=f'{self.outdir}/{samp}_gene_quant'
        cmd=[utils.selectENV('VC-Quantify')]
        cmd.extend(
            ['salmon quant --validateMappings','-i',salmon_idx,
            '-l A','-p',self.threads,'--meta',
            '-1',self.fastqs[0],'-2',self.fastqs[1],'-o',wkdir,'\n']
        )
        return cmd
