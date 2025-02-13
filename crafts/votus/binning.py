from ..general import utils
from ..data.bioseq import VirSeq

class VirMAG(VirSeq):
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    
    def filtCtg(self, checkv=None):
        cmd = []
        ckeckv_dir = ''
        if not checkv:
            cmd, ckeckv_dir = self.checkv(self.fasta)
        return cmd

    def vrhyme(self, gene_fa=None, protein_fa=None, coverage_tsv=None):
        cmd = [utils.selectENV('VC-vRhyme')]
        wkdir = f'{self.wkfile_dir}/vrhyme'
        cmd.extend(
            ['vRhyme', '-i', self.fasta, '-g', gene_fa, '-p', protein_fa, '-c', coverage_tsv, '-t', self.threads, '-o', wkdir, '\n']
        )
        return cmd, wkdir

    def Binning(self, gene_fa: str, protein_fa: str, coverage_tsv: str):
        cmd = []
        if not (gene_fa or protein_fa):
            cmd, protein_fa = self.genePred()
            gene_fa = protein_fa.replace('.faa', '.ffn')
        if not coverage_tsv:
            pass
        cmd = vrhyme(gene_fa, protein_fa, coverage_tsv)
        shell = f'{self.shell_dir}/{self.name}_viral_binning.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results

