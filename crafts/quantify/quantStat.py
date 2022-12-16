from .multiQuant import multiCount
from ..general import cmdExec,general

class AbdStat(multiCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def mergeAbd(self): #Heatmap for contigs abundance
        abd=f'{self.outdir}/all_merged.tpm'
        cmd=['merge_tpms.pl',self.samp_info,self.outdir,'tpm\n',
            'pheatmap_for_abd.R',abd,self.samp_info,self.outdir,'tpm\n']
        return cmd,abd
    def sizeAbdPlot(self,abd:str):
        tmp_cmd,fasta_stat=self.sizeGC() #calculate the size of each contig
        cmd=tmp_cmd
        abd=f'{self.outdir}/all_merged'
        sum_abd=f'{self.outdir}/all_sum_abd.xls'
        sum_len_abd=f'{self.outdir}/all_sum_len_abd.xls'
        cmd.extend(
            ['fasta_size_gc.py',self.fasta,'>',fasta_stat,'\n',
            'sum_abd_by_seq.py',abd,self.outdir,'\n',
            'linkTab.py',fasta_stat,sum_abd,'left Contig',sum_len_abd,'\n',
            'fa_length_tpm_scatter.R',sum_len_abd,self.outdir,'\n']
        )
        return cmd
    def taxaAbd(self,abd:str,taxa_anno:str):
        modi_tax_anno=f'{self.outdir}/DemoVir_assignments.txt'
        cmd=["sed '1s/Sequence_ID/Contig/'",tax_anno,'>',modi_tax_anno,'\n',
            'linkTab.py',abd,'left Contig',modi_tax_anno,'\n',
            'sum_abd_by_taxa.py',abd,self.outdir,'\n',
            'barplot_for_taxa_tpm.R',tax_tpm,self.outdir,'\n']
        return cmd
    def diversity(self,abd:str):        
        cmd=['alpha_diversity.R',abd,self.samp_info,self.outdir,'\n',
            'NMDS.R',abd,self.samp_info,self.outdir,'\n']
        return cmd
    def QuantStat(self,taxa_anno:str):
        self.countBySamp()
        cmd=[self.envs]
        tmp_cmd,abd=self.mergeAbd()
        cmd.extend(tmp_cmd)
        cmd.extend(self.sizeAbdPlot(abd))
        if taxa_anno: cmd.extend(self.taxaAbd(abd,self.taxa_anno))
        cmd.extend(self.diversity(abd))
        shell=f'{self.outdir}/{self.name}_count.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
