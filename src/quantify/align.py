#!/usr/bin/env python3

import sys,os
from ..general import cmdExec,general
from ..fastqc.readsQC import Reads

class VirAln(Reads):
    '''
    '''
    def __init__(self,fq1='',fq2='',fasta='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
    def mkIdx(self):
        "Make bwa index for votus."
        cmd=[self.envs]
        wkdir=f'{self.wkdir}/0.index'
        general.mkdir(wkdir)
        votus_lnk=f'{wkdir}/all_votus.fa'
        cmd.extend(['ln -s',self.votus,votus_lnk,'\n'])
        shell=f'{self.wkdir}/bwa_index.sh'
        cmd.extend(['bwa index -a bwtsw',votus_lnk,'-p',self.bwa_idx,'\n'])
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def alnReads(self,samp:str,bwa_idx:str,wkdir:str):
        '''
        Align the reads to the vOTUs for each sample.
        '''
        cmd=[self.envs]
        raw_bam=f'{wkdir}/{samp}.raw.bam'
        sort_bam=f'{wkdir}/{samp}.sort.bam'
        fastqs=self.sampDict[samp][2].split(',')
        cmd.extend(
            ['bwa mem','-t 28',bwa_idx,fastqs[0],fastqs[1],
             '|samtools view','-o',raw_bam,'-@ 28 -b -S\n',
             'samtools sort',raw_bam,'-o',sort_bam,'-@ 28\n']
        )
        shell=f'{outdir}/0.bwa/{samp}_bwa.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def AlnBySamp(self):
        wkdir=f'{self.wkdir}/1.bwa'
        general.mkdir(wkdir)
        results=mkIdx()
        for samp in self.sampDict.keys():
            results+=self.alnReads(samp,bwa_idx,wkdir)
        return results
