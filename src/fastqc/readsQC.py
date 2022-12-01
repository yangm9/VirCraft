import os
from ..general import cmdExec,general
from ..config.cfgInfo import VirCfg

class Reads:
    '''
    FastQ processing class.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,fq1='',fq2='',outdir=''):
        self.fastqs=[fq1,fq2]
        self.basename_fq1=os.path.basename(self.fastqs[0])
        self.basename_fq2=os.path.basename(self.fastqs[1])
        self.outdir=os.path.abspath(outdir)
        self.samp=self.basename_fq1.replace('_1.fastq','')
        self.samp=self.samp.replace('_1.fq','')
        general.mkdir(self.outdir)

class qualCtrl(Reads):
    def __init__(self,fq1='',fq2='',outdir=''):
        super().__init__(fq1,fq2,outdir)
    def filtReads(self):
        wkdir=f'{self.outdir}/fastp'
        out_fq1=f'{wkdir}/{self.basename_fq1}'
        out_fq2=f'{wkdir}/{self.basename_fq2}'
        filt_rpts=f'{out_fq1}_report.html'
        out_fq_list=f'{out_fq1}_list.txt'
        general.mkdir(wkdir)
        cmd=['fastp','-i',fq1,'-o',out_fq1,'-I',fq2,'-O',out_fq2,
            '-w 16 -q 20 -u 20 -g -c -W 5 -3 -l 50','-h',filt_rptsi,'\n',
            'ls',out_fq1,out_fq2,'>',out_fq_list,'\n']
        return cmd
    def fastUniq(self):
        wkdir=f'{self.outdir}/fastuniq'
        in_fq_list=f'{self.outdir}/fastp/{self.basename_fq1}_list.txt'
        out_fq1=f'{wkdir}/{self.basename_fq1}'
        out_fq2=f'{wkdir}/{self.basename_fq2}'
        general.mkdir(wkdir)
        cmd=['fastp','-i',out_fq_list,'-o',out_fq1,'-p',out_fq2,'\n']
        return cmd
    def decontaminate(self,fq1,fq2):
        wkdir=f'{self.outdir}/decontaminate'



        
        
        
        
    

