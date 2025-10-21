from ..general import utils
from ..identify.viralDetectors import VirDetectTools

class VirMAG(VirDetectTools):
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    
    def filt_ctgs_for_binning(self, checkv=None):
        cmd = [utils.selectENV('VC-General')]
        if not checkv:
            outdir = self.outdir
            self.outdir = self.wkfile_dir
            cmd, checkv_dir = self.checkv(self.fasta)
            checkv = f'{checkv_dir}/quality_summary.tsv'
            self.outdir = outdir
        vctg_for_binning = f'{self.wkfile_dir}/{self.name}_vctg_for_binning.fa'
        cmd.extend(
            ['keep_vctg_for_binning.pl', checkv, self.fasta, vctg_for_binning, '\n']
        )
        return cmd, vctg_for_binning

    def vrhyme(self, fasta_for_binning: str, coverage_tsv=None, fastqs=None):
        cmd = [utils.selectENV('VC-vRhyme')]
        wkdir = f'{self.wkfile_dir}/vrhyme'
        if coverage_tsv:
            sub_cmd =  f'-c {coverage_tsv}'
        elif fastqs:
            input_arg = f'-r {fastqs}'
        else:
            raise ValueError('Error: Either coverage_tsv or fastqs must be provided')
        cmd.extend(
            ['vRhyme', '-i', self.fasta, '-g', gene_fa, '-p', protein_fa, sub_cmd, '-t', self.threads, '-o', wkdir, '\n']
        )
        return cmd, wkdir

    def vMAGFilt(self):
        pass

    def Binning(self, cov_input=None, checkv=None, unrun=False):
        cmd, vctg_for_binning = self.filtCtg(checkv)
        tmp_cmd, gene_fa, protein_fa = self.genePred()
        cmd.extend(tmp_cmd)
        tmp_cmd, vrhyme_dir = self.vrhyme(vctg_for_binning, coverage_tsv)
        cmd.extend(tmp_cmd)
        tmp_cmd, = vMAGFilt()
        cmd.extend(tmp_cmd)
        shell = f'{self.shell_dir}/{self.name}_vctg_binning.sh'
        utils.printSH(shell, cmd)
        rcode = 0 if unrun else utils.execute(shell)
        return rcode
