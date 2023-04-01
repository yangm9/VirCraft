from ..general import utils
from .multiQuant import multiGeneCount

class GeneAbdStat(multiGeneCount):
    '''
    Gene abundance main class
    '''
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
        self.parallel=str(parallel)
    def mergeAbd(self):
        abd=f'{self.outdir}/all_merged_gene.sf'
        cmd=['merge_tpms.pl',self.samp_info,self.outdir,'tpm Gene\n']
        return cmd,abd
    def QuantStat(self,batch_size):
        self.geneCountBySamp()
        cmd=[self.envs]
        cmd.extend(
            ['multithreads.pl',self.outdir,'gene_count.sh',batch_size,'\n']
        )
        tmp_cmd,abd=self.mergeAbd()
        cmd.extend(tmp_cmd)
        cmd.append('rm -rf *_gene_count.sh*\n')
        shell=f'{self.outdir}/{self.name}_gene_quant.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(cmd)
        return results
