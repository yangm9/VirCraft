import os
from ..general import utils
from ..data.bioseq import Seq
from ..data.bioseq import CDS
from .alnQuant import VirCount
from .alnQuant import GeneCount

#Running virus or gene abundance analysis in batch

class multiVirCount(Seq):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.samp_info=os.path.abspath(samp_info)
        self.groups,self.sampDict=self.readSampInfo(self.samp_info)
        self.threads=str(threads)
    def virCountBySamp(self):
        idx_cmd,bwa_idx=self.mkBwaIdx()
        for samp in self.sampDict.keys():
            cmd=[utils.selectENV('VC-General')]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=VirCount(fq1,fq2,self.wkdir,self.threads)
            cmd.extend(Count.bwa(samp,bwa_idx))
            cmd.extend(Count.coverm(samp))
            shell=f'{self.shelldir}/{samp}_viral_count.sh'
            utils.printSH(shell,cmd)
        return idx_cmd

class multiGeneCount(CDS):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
        self.samp_info=os.path.abspath(samp_info)
        self.groups,self.sampDict=self.readSampInfo(self.samp_info)
        self.threads=str(threads)
    def geneCountBySamp(self):
        cmd=[utils.selectENV('VC-General')]
        idx_cmd,salmon_idx=self.mkSalmonIdx
        for samp in self.sampDict.keys():
            cmd=[self.envs]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=GeneCount(fq1,fq2,self.outdir,self.threads)
            cmd.extend(Count.salmon(samp,salmon_idx))
            shell=f'{self.outdir}/{samp}_gene_count.sh'
            utils.printSH(shell,cmd)
        return idx_cmd 
