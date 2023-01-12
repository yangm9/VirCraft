
from ..general import cmdExec,general
from ..config.config import Seq

class AMGs(Seq):
    '''
    Call AMGs.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def dramv(self):
        vs2dir=f'{self.outdir}/vs2'
        wkdir=f'{self.outdir}/dramv'
        vs2_dramv_fa=f'{vs2dir}/for-dramv/final-viral-combined-for-dramv.fa'
        vs2_dramv_tab=f'{vs2dir}/for-dramv/viral-affi-contigs-for-dramv.tab'
        anno_tsv=f'{wkdir}/annotations.tsv'
        distill_dir=f'{wkdir}/distilled'
        cmd=['virsorter run --prep-for-dramv','-w',wkdir,
            '-i',self.fasta,'-j',self.threads,'all\n',
            'DRAM-v.py annotate','--threads',self.threads,
            '-i',vs2_dramv_fa,'-v',vs2_dramv_tab,'-o',wkdir,'\n',
            'DRAM-v.py distill','-i',anno_tsv,'-o',distill_dir,'\n']
        return cmd
    def vibrant(self):
        wkdir=f'{self.outdir}/vibrant'
        return cmd
