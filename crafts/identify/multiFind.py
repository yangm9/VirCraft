import os
from ..general import utils
from ..identify.viridsop import VirScan

class MultiTools(VirScan):
    '''
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def deepvirfinder(self,cutoff:str):
        cmd=[utils.selectENV('deepvirfinder')]
        wkdir=f'{self.outdir}/deepvirfinder'
        utils.mkdir(wkdir)
        dvf=self.confDict['DeepVirFinder']
        models=os.path.dirname(dvf)+'/models'
        cmd.extend(
            ['python',dvf,'-i',self.fasta,'-m',models,'-o',wkdir,
            '-c',self.threads,'-l',cutoff,'\n']
        )
        return cmd,wkdir
    def vibrant(self):
        cmd=[utils.selectENV('vibrant')]
        wkdir=f'{self.outdir}/VIBRANT_{self.name}'
        cmd.extend(
            ['VIBRANT_run.py','-i',self.fasta,
            '-t',self.threads,'-folder',self.outdir,'\n']
        )
        return cmd,wkdir
