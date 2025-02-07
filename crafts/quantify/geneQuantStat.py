from ..general import utils
from .multiQuant import multiGeneCount

class GeneAbdStat(multiGeneCount):
    '''
    Gene abundance main class.
    '''
    def __init__(self, samp_info=None, fasta=None, outdir=None, threads=8):
        "The FastA file should be coding sequence file (*.ffn)."
        super().__init__(samp_info, fasta, outdir, threads)
    def mergeAbd(self):
        abd = f'{self.wkdir}/all_merged_gene.sf'
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(['paste_abundance_files.py', self.samp_info, self.wkdir, 'sf Gene\n'])
        return cmd, abd
    def QuantStat(self, unrun=False):#,batch_size):
        cmd = self.geneCountBySamp()
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(
            ['multithreads.pl', self.shelldir, 'gene_count.sh', str(self.BATCH_SIZE), '\n']
        )
        tmp_cmd,abd = self.mergeAbd()
        cmd.extend(tmp_cmd)
        #cmd.append('rm -rf *_gene_count.sh*\n')
        shell = f'{self.shelldir}/{self.name}_gene_quant.sh'
        results = ''
        utils.printSH(shell, cmd)
        if not unrun: results = utils.execute(shell)
        return results
