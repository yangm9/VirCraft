import os
from ..general import utils
from ..data.bioseq import Seq
from .alnQuant import VirCount

class multiVirCount(Seq):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.samp_info=os.path.abspath(samp_info)
        self.groups,self.sampDict=self.readSampInfo(self.samp_info)
        self.threads=str(threads)
    def virCountBySamp(self):
        bwa_idx,_=self.mkBwaIdx
        results=''
        for samp in self.sampDict.keys():
            cmd=[self.envs]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=VirCount(fq1,fq2,self.outdir,self.threads)
            cmd.extend(Count.bwa(samp,bwa_idx))
            cmd.extend(Count.coverm(samp))
            shell=f'{self.outdir}/{samp}_viral_count.sh'
            utils.printSH(shell,cmd)
            results+=utils.execute(cmd)
        return results

class multiGeneCount(multiVirCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def geneCountBySamp(self):
        salmon_idx,_=self.mkSalmonIdx
        results=''
        for samp in self.sampDict.keys():
            cmd=[self.envs]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=GeneCount(fq1,fq2,self.outdir,self.threads)
            cmd.extend(Count.salmon(samp,salmon_idx))
            shell=f'{self.outdir}/{samp}_gene_count.sh'
            utils.printSH(shell,cmd)
            results+=utils.execute(cmd)
        return results
