import os
from ..general import utils
from ..identify.viridsop import VirScan

class MultiTools(VirScan):
    '''
    Generate command for deepvirfinder and VIBRANT.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def deepvirfinder(self,cutoff:str):
        cmd=[utils.selectENV('VC-DeepVirFinder')]
        wkdir=f'{self.outdir}/deepvirfinder'
        utils.mkdir(wkdir)
        cmd.extend(
            ['dvf.py','-i',self.fasta,'-m',self.confDict['DeepVirFinderDB'],
            '-o',wkdir,'-c',self.threads,'-l',cutoff,'\n']
        )
        return cmd,wkdir
    def vibrant(self,cutoff:str):
        cmd=[utils.selectENV('VC-VIBRANT')]
        wkdir=f'{self.outdir}/VIBRANT_{self.name}'
        VIBRANT_DB_databases=self.confDict['VIBRANTDB']+'/databases'
        VIBRANT_DB_files=self.confDict['VIBRANTDB']+'/files'
        cmd.extend(
            ['VIBRANT_run.py','-i',self.fasta,'-l',cutoff,
            '-t',self.threads,'-folder',self.outdir,
            '-d',VIBRANT_DB_databases,'-m',VIBRANT_DB_files,'\n']
        )
        return cmd,wkdir
