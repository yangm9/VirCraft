from ..general import utils
from .viralDetectors import VirDetectTools

class VirMAG(VirDetectTools):
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    
    def vrhyme(self, coverage_mode: str, coverage_files: str): # mode: 'pe_fastq' | 'se_fastq' | 'sam' | 'bam' | 'coverage'
        cmd = [utils.selectENV('VC-vRhyme')]
        wkdir = f'{self.wkfile_dir}/vrhyme'
        
        tmp_cmd, gene_fa, protein_fa = self.genePred()
        cmd.extend(tmp_cmd)
        cov_file_txt = ' '.join(coverage_files)
        sub_cmd = ''
        if coverage_mode == 'pe_fastq':
            sub_cmd = f'-r {cov_file_txt}'
        elif coverage_mode == 'se_fastq':
            sub_cmd = f'-u {cov_file_txt}'
        elif coverage_mode == 'sam':
            sub_cmd = f'-s {cov_file_txt}'
        elif coverage_mode == 'bam':
            sub_cmd = f'-b {cov_file_txt}'
        elif coverage_mode == 'cov':
            sub_cmd = f'-c {cov_file_txt}'
        else:
            raise ValueError('Error: mode must be "pe_fastq", "se_fastq", "sam", "bam" or "cov".')
        
        cmd.extend(
            ['vRhyme', '-i', self.fasta, '-g', gene_fa, '-p', protein_fa, sub_cmd, '-t', self.threads, '-o', wkdir, self.confDict['vRhymeOpts'], '\n']
        )
        return cmd, wkdir

    def vMAGFilt(self):
        cmd = [utils.selectENV('VC-CAT')]
        return cmd

    def vMAGQC(self):
        return 0
    def Binning(self, coverage_mode=None, coverage_files=None, unrun=False):
        cmd, vrhyme_dir = self.vrhyme(coverage_mode, coverage_files)
        tmp_cmd = self.vMAGFilt()
        cmd.extend(tmp_cmd)
        shell = f'{self.shell_dir}/{self.name}_vctg_binning.sh'
        utils.printSH(shell, cmd)
        rcode = 0 if unrun else utils.execute(shell)
        return rcode
