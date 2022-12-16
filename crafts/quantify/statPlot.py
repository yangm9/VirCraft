import os
from .multiQuant import multiCount

class AbdStat(multiCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def calcAbd(self): #Heatmap for contigs abundance
        abd=f'{self.outdir}/all_merged.tpm'
        cmd.extend(
            ['merge_tpms.pl',self.samp_info,self.outdir,'tpm\n',
            'pheatmap_for_abd.R',abd,self.samp_info,self.outdir,'tpm\n']
        )
        return cmd,abd
    def sizeAbdPlot(self,abd:str):
        tmp_cmd,fasta_stat=self.sizeGC()
        cmd.extend(tmp_cmd) #calculate the size of each contig
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

