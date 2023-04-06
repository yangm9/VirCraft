import os
from ..general import utils
from ..identify.viridsop import VirScan

class MultiTools(VirScan):
    '''
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def deepvirfider(self,cutoff:str):
        cmd=[utils.selectENV('deepvirfider')]
        wkdir=f'{self.outdir}/deepvirfinder'
        cmd.extend(
            ['python',dvf,'-i',self.fasta,'-m',models,'-o',wkdir,
            '-c',self.threads,'-l',cutoff,'\n']
        )
        return cmd
    def vibrant(self):
        cmd=[utils.selectENV('vibrant')]
        wkdir=f'{self.outdir}/vibrant'
        cmd=['VIBRANT_run.py','-i',self.fasta,
            '-t',self.threads,'-folder',wkdir,'\n']
        return cmd
    def vFilter(self):
        
    def Identify(self):
        cmd.extend(tmp_cmd)
        #Step 4 DRAMv annotation
        tmp_cmd=self.annotate(wkdir) 
        cmd.extend(tmp_cmd)
        #Step 5 Curation
        tmp_cmd=self.curate()
        cmd.extend(tmp_cmd)
        #Generate shell and exeute it
        shell=f'{self.outdir}/{self.name}_find_vir.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(cmd)
        return results
