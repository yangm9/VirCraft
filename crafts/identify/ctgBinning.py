from ..general import utils
from .viralDetectors import VirDetectTools

class VirMAG(VirDetectTools):
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    
    def vrhyme(self, vir_ctgs_fa: str, mode: str, input_files: str): # mode: 'pe_fastq' | 'se_fastq' | 'sam' | 'bam' | 'coverage'
        cmd = [utils.selectENV('VC-vRhyme')]
        wkdir = f'{self.wkfile_dir}/vrhyme'
        
        tmp_cmd, gene_fa, protein_fa = self.genePred()
        cmd.extend(tmp_cmd)
        
        sub_cmd = ''
        if mode == 'pe_fastq':
            sub_cmd = ['-r', input_files]
        elif se_fastqs:
            sub_cmd = ['-u', input_files]
        elif sams:
            sub_cmd = ['-s', input_files]
        elif bams:
            sub_cmd = ['-b', input_files]
        elif coverage_tsv:
            sub_cmd = ['-c', input_files]
        else:
            raise ValueError('Error: mode must be "pe_fastqs", "se_fastqs", "sams", "bams" or "coverage_tsv".')
        
        cmd.extend(
            ['vRhyme', '-i', self.fasta, '-g', gene_fa, '-p', protein_fa, sub_cmd, '-t', self.threads, '-o', wkdir, '\n']
        )
        return cmd, wkdir

    def vMAGFilt(self):
        cmd = [utils.selectENV('VC-CAT')]
        pass

    def Binning(self, mode=None, fsbc_files=None, checkv=None, unrun=False):
        cmd, vrhyme_dir = self.vrhyme(vctg_for_binning, mode, fsbc_files)
        cmd.extend(tmp_cmd)
        tmp_cmd = vMAGFilt()
        cmd.extend(tmp_cmd)
        shell = f'{self.shell_dir}/{self.name}_vctg_binning.sh'
        utils.printSH(shell, cmd)
        rcode = 0 if unrun else utils.execute(shell)
        return rcode
