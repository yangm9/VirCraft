import os
import re
import sys
from ..general import general
from ..general import cmdExec

class Reads(VirCfg):
    '''
    FastQ processing class.
    '''
    envs=general.selectENV('VirCraft')
    postfixes=[
        '_1.fastq','_1.fastq.gz','_1.fq','_1.fq.gz',
        '_R1.fastq','_R1.fastq.gz','_R1.fq','_R1.fq.gz',
        '.R1.fastq','.R1.fastq.gz','.R1.fq','.R1.fq.gz'
    ]
    def __init__(self,fq1='',fq2='',outdir='',*args,**kwargs):
        super().__init__()
        fq1=os.path.abspath(fq1)
        fq2=os.path.abspath(fq2)
        self.fastqs=[fq1,fq2]
        self.basename_fq1=os.path.basename(self.fastqs[0])
        self.basename_fq2=os.path.basename(self.fastqs[1])
        self.outdir=os.path.abspath(outdir)
        self.samp=self.getSampName
        general.mkdir(self.outdir)
    @property
    def getSampName(self):
        samp=''
        for post in self.postfixes:
            if self.basename_fq1.endswith(post):
                samp=self.basename_fq1.replace(post,'')
        return samp

class Seq(VirCfg):
    '''
    Fasta processing class.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,fasta='',outdir='',*args,**kwargs):
        super().__init__()
        basename_fa=os.path.basename(fasta)
        self.name=os.path.splitext(basename_fa)[0]
        self.fasta=os.path.abspath(fasta)
        self.outdir=os.path.abspath(outdir)
        general.mkdir(self.outdir)
    @property
    def mkBwaIdx(self):
        "Make bwa index for votus."
        cmd=[self.envs]
        idx=f'{self.outdir}/{self.name}BWAIDX'
        cmd.extend(['bwa index -a bwtsw',self.fasta,'-p',idx,'\n'])
        shell=f'{self.outdir}/{self.name}_bwaidx.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return idx,results
    def sizeGC(self): #Plot scatter for contigs size (x) and GC content (y)
        cmd=[self.envs]
        fasta_stat=f'{self.outdir}/fasta_stat.xls'
        cmd.extend(
            ['fasta_size_gc.py',self.fasta,'>',fasta_stat,'\n',
            'variables_scatter.R',fasta_stat,'Length~GC',self.outdir,'\n']
        )
        return cmd,fasta_stat
    def genePred(self):
        cmd=[self.envs]
        wkdir=f'{self.outdir}/prodigal'
        general.mkdir(wkdir)
        temp_orfs_faa=f'{wkdir}/temp.orfs.faa'
        orfs_ffn=f'{wkdir}/{self.name}_votus.ffn'
        orfs_faa=f'{wkdir}/{self.name}_votus.faa'
        temp_orfs_ffn=f'{wkdir}/temp.orfs.ffn'
        temp_txt=f'{wkdir}/temp.txt'
        cmd=['prodigal','-i',self.fasta,'-a',temp_orfs_faa,
            '-d',temp_orfs_ffn,'-m','-o',temp_txt,'-p meta -q\n',
            'cut -f 1 -d \" \"',temp_orfs_faa,'>',orfs_faa,'\n',
            'cut -f 1 -d \" \"',temp_orfs_ffn,'>',orfs_ffn,'\n',
            f'rm -f {wkdir}/temp.*\n']
        return cmd,orfs_faa

class VirSeq(Seq):
    def __init__(self,fasta='',outdir='',*args,**kwargs):
        super().__init__(fasta,outdir,*args,**kwargs)
    def checkv(self,in_fa:str):
        cmd=[self.envs]
        wkdir=f'{self.outdir}/checkv'
        general.mkdir(wkdir)
        cmd=['checkv','end_to_end',in_fa,wkdir,
            '-d',self.confDict['CheckvDB'],
            '-t',self.threads,'\n']
        provir_fna=f'{wkdir}/proviruses.fna'
        vir_fna=f'{wkdir}/viruses.fna'
        out_fa=f'{wkdir}/combined.fna'
        cmd.extend(['cat',provir_fna,vir_fna,'>',out_fa,'\n'])
        return cmd,out_fa

class ORF(Seq):
    def __init__(self,orfs='',outdir='',*args,**kwargs):
        super().__init__(orfs,outdir,*args,**kwargs)
        self.orfs=self.fasta
    @property
    def mkSalmonIdx(self):
        cmd=[self.envs]
        idx=f'{self.outdir}/SalmonIdx' # A directory
        general.mkdir(wkdir)
        cmd.extend(
            ['salmon index','-p 9 -k 31','-t',self.orfs,'-i',wkdir,'\n']
        )
        shell=f'{self.outdir}/{self.name}_salmonidx.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return idx,results
