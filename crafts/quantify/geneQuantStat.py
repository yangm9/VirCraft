from .multiQuant import multiGeneCount
from ..general import cmdExec,general

class GeneAbdStat(multiGeneCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def mergeAbd(self):
        abd=f'{self.outdir}/all_merged_gene.tpm'
        cmd=['merge_tpms.pl',self.samp_info,self.outdir,'tpm\n']
        return cmd,abd
