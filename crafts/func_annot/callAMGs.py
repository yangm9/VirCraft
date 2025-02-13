from ..general import utils
from ..identify.viralDetectors import VirDetectTools
from .geneAnnot import GeneFunc

class AMGs(VirDetectTools):
    '''
    Call AMGs.
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir, threads)
        self.threads = str(threads)
    def dramv(self, vs2_dramv_fa, vs2_dramv_tab):
        cmd = [utils.selectENV('VC-DRAMv')]
        wkdir = f'{self.wkfile_dir}/dramv'
        anno_tsv = f'{wkdir}/annotations.tsv'
        distill_dir = f'{wkdir}/distilled'
        cmd.extend(
            ['rm -rf',wkdir,'\n',
            'DRAM-v.py annotate', '--threads', self.threads, '-i', vs2_dramv_fa, '-v', vs2_dramv_tab, '-o', wkdir, '\n',
            'DRAM-v.py distill', '-i', anno_tsv, '-o', distill_dir, '\n']
        )
        return cmd, wkdir
    def mergeResults(self, dmvdir: str, vbdir: str):
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(['vb_dmv_annot_merge.py', dmvdir, vbdir, self.wkfile_dir, '\n'])
        return cmd
    def annotBatchSH(self, orf_faa: str):
        outdir = self.outdir
        self.outdir = self.wkfile_dir
        cmd = self.eggnogAnno(orf_faa) # invoke eggnogAnno() method from GeneFunc
        shell = f'{self.shell_dir}/{self.name}_eggnog_anno.sh'
        utils.printSH(shell, cmd)
        cmd = self.keggAnno(orf_faa) # invoke keggAnno() method from GeneFunc
        self.outdir = outdir
        shell = f'{self.shell_dir}/{self.name}_kegg_anno.sh'
        utils.printSH(shell, cmd)
        return 0
    def amgBatchSH(self, vs2_dramv_fa: str, vs2_dramv_tab: str):
        cmd, dmvdir = self.dramv(vs2_dramv_fa, vs2_dramv_tab)
        shell = f'{self.shell_dir}/{self.name}_dramv_amg.sh'
        utils.printSH(shell, cmd)
        outdir = self.outdir
        self.outdir = self.wkfile_dir
        cmd, vbdir = self.vibrant()
        self.outdir = outdir #self.outdir need to be changed back to its original value
        shell = f'{self.shell_dir}/{self.name}_vibrant_amg.sh'
        utils.printSH(shell, cmd)
        return dmvdir, vbdir
    def annotAMGs(self, unrun=False):
        #step1 VS2
        cmd, vs2dir = self.virsorter(self.fasta, 1, min_score=0.5, min_length=1500)
        #step2 for DRAM-v and VIBRANT 
        vs2_dramv_fa = f'{vs2dir}/for-dramv/final-viral-combined-for-dramv.fa'
        vs2_dramv_tab = f'{vs2dir}/for-dramv/viral-affi-contigs-for-dramv.tab'
        dmvdir, vbdir = self.amgBatchSH(vs2_dramv_fa, vs2_dramv_tab)
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(['multithreads.pl', self.shell_dir, 'amg.sh 2\n'])
        #step3 prodigal
        self.fasta = vs2_dramv_fa
        tmp_cmd, orf_faa = self.genePred()
        cmd.extend(tmp_cmd)
        #step4 eggnog and kegg annotation
        self.annotBatchSH(orf_faa)
        cmd.extend(['multithreads.pl', self.shell_dir, 'anno.sh 2\n'])
        #step5 merge all AMGs results
        cmd.extend(self.mergeResults(dmvdir, vbdir))
        shell = f'{self.shell_dir}/{self.name}_call_amgs.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
