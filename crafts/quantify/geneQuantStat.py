from .multiQuant import multiGeneCount

class GeneAbdStat(multiGeneCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def mergeAbd(self):
        abd=f'{self.outdir}/all_merged_gene.sf'
        cmd=['merge_tpms.pl',self.samp_info,self.outdir,'tpm\n']
        return cmd,abd
    def QuantStat(self):
        self.geneCountBySamp()
        cmd=[self.envs]
        tmp_cmd,abd=self.mergeAbd()
        cmd.extend(tmp_cmd)
        shell=f'{self.outdir}/{self.name}_gene_count.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(cmd)
        return results
