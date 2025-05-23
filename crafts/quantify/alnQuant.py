from ..general import utils
from ..data.bioseq import Reads

class VirCount(Reads):
    def __init__(self, fq1=None, fq2=None, outdir=None, threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads = str(int(threads) // self.BATCH_SIZE)
    def bwa(self, samp: str, bwa_idx: str):
        '''
        Align the reads to the vOTUs for each sample.
        '''
        raw_bam = f'{self.wkfile_dir}/{samp}.raw.bam'
        sort_bam = f'{self.wkfile_dir}/{samp}.sort.bam'
        cmd = [utils.selectENV('VC-Quantify')]
        cmd.extend(
            ['bwa mem', '-t', self.threads, bwa_idx, self.fastqs[0], self.fastqs[1], '|samtools view', '-o', raw_bam, '-@ 28 -b -S\n',
            'samtools sort', raw_bam, '-o', sort_bam,'-@ 28\n', 'rm -f', raw_bam, '\n', 'samtools index', sort_bam, '\n']
        )
        return cmd
    def coverm(self, samp: str, method='mean'):
        cov = f'{self.wkfile_dir}/{samp}.cov'
        sort_bam = f'{self.wkfile_dir}/{samp}.sort.bam'
        cmd = [utils.selectENV('VC-Quantify')]
        tmp_cmd = ['coverm contig', '-m', method, '-b', sort_bam, '-t', self.threads, self.confDict['CoverMOpts'], '>', cov, '\n']
        if method == 'metabat':
            tmp_cmd.remove(self.confDict['CoverMOpts'])
        cmd.extend(tmp_cmd)
        return cmd

class GeneCount(Reads):
    def __init__(self, fq1=None, fq2=None, outdir=None, threads=8):
        super().__init__(fq1, fq2, outdir)
        self.threads = str(int(threads) // self.BATCH_SIZE)
    def salmon(self, samp: str, salmon_idx: str):
        wkdir = f'{self.wkfile_dir}/{samp}_gene_quant'
        cmd = [utils.selectENV('VC-Quantify')]
        cmd.extend(
            ['salmon quant --validateMappings','-i', salmon_idx, '-l A', '-p', self.threads, '--meta', '-1', self.fastqs[0], '-2', self.fastqs[1], '-o', wkdir, '\n']
        )
        return cmd
