from ..general import utils
from ..identify.viridsop import VirScan
from ..identify.multiFind import MultiTools
from .geneAnnot import GeneFunc

class AMGs(VirScan):
    '''
    Call AMGs.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
        self.threads=str(threads)
    def dramv(self,vs2_dramv_fa,vs2_dramv_tab):
        cmd=[utils.selectENV('viral-id-sop')]
        wkdir=f'{self.outdir}/dramv'
        anno_tsv=f'{wkdir}/annotations.tsv'
        distill_dir=f'{wkdir}/distilled'
        cmd.extend(
            ['DRAM-v.py annotate','--threads',self.threads,
            '-i',vs2_dramv_fa,'-v',vs2_dramv_tab,'-o',wkdir,'\n',
            'DRAM-v.py distill','-i',anno_tsv,'-o',distill_dir,'\n']
        )
        return cmd,wkdir
    def mergeResults(self,dmvdir,vbdir):
        cmd=[utils.selectENV('VirCraft')]
        cmd.extend(['AMGs_filter.py',dmvdir,vbdir,self.outdir,'\n'])
        return cmd
    def annotBatchSH(self,orf_faa):
        GF=GeneFunc(
            fasta=self.fasta,
            outdir=self.outdir,
            threads=self.threads
        )
        cmd=[self.envs]
        cmd.extend(GF.eggnogAnno(orf_faa))
        shell=f'{self.outdir}/{self.name}_eggnog_anno.sh'
        utils.printSH(shell,cmd)
        cmd=[utils.selectENV('kofamscan')]
        cmd.extend(GF.keggAnno(orf_faa))
        shell=f'{self.outdir}/{self.name}_kegg_anno.sh'
        utils.printSH(shell,cmd)
        return 0
    def amgBatchSH(self,vs2_dramv_fa,vs2_dramv_tab):
        cmd,dmvdir=self.dramv(vs2_dramv_fa,vs2_dramv_tab)
        shell=f'{self.outdir}/{self.name}_dramv_amg.sh'
        utils.printSH(shell,cmd)
        MT=MultiTools(
            fasta=vs2_dramv_fa,
            outdir=self.outdir,
            threads=self.threads
        )
        cmd,vbdir=MT.vibrant()
        shell=f'{self.outdir}/{self.name}_vibrant_amg.sh'
        utils.printSH(shell,cmd)
        return dmvdir,vbdir
        
    def annotAMGs(self):
        #step1 VS2
        cmd=[utils.selectENV('viral-id-sop')]
        tmp_cmd,vs2dir=self.virsorter(self.fasta,1)
        cmd.extend(tmp_cmd)
        #step2 for DRAM-v and VIBRANT 
        vs2_dramv_fa=f'{vs2dir}/for-dramv/final-viral-combined-for-dramv.fa'
        vs2_dramv_tab=f'{vs2dir}/for-dramv/viral-affi-contigs-for-dramv.tab'
        dmvdir,vbdir=self.amgBatchSH(vs2_dramv_fa,vs2_dramv_tab)
        cmd.extend(['multithreads.pl',self.outdir,'amg.sh 2\n'])
        #step3 prodigal
        self.fasta=vs2_dramv_fa
        tmp_cmd,orf_faa=self.genePred()
        cmd.extend(tmp_cmd)
        #step4 eggnog and kegg annotation
        self.annotBatchSH(orf_faa)
        cmd.extend(['multithreads.pl',self.outdir,'anno.sh 2\n'])
        #step5 merge all AMGs results
        cmd.extend(self.mergeResults(dmvdir,vbdir))
        shell=f'{self.outdir}/{self.name}_call_amgs.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(cmd)
        return results
