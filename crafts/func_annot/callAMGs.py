
from ..general import cmdExec,general
from ..config.config import Seq

class AMGs(Seq):
    '''
    Call AMGs.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
