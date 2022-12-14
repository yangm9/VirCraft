import os
from ..general import cmdExec,general
from .alnQuant import VirCount
from ..config.config import Seq

class multiCount(Seq):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.samp_info=os.path.abspath(samp_info)
        self.groups,self.sampDict=self.readSampInfo(self.samp_info)
        self.threads=str(threads)
    def countBySamp(self):
        bwa_idx,_=self.mkBwaIdx
        results=''
        for samp in self.sampDict.keys():
            cmd=[self.envs]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=VirCount(fq1,fq2,self.outdir,self.threads)
            cmd.extend(Count.aln(samp,bwa_idx))
            cmd.extend(Count.coverm(samp))
            shell=f'{self.outdir}/{samp}_count.sh'
            general.printSH(shell,cmd)
            results+=cmdExec.execute(cmd)
        return results
    def statPlot(self,tax_anno=''):
        cmd=[self.envs]
        cmd.extend(['merge_tpms.pl',self.samp_info,self.outdir,'tpm\n'])
        modi_tax_anno=general.insLable(tax_anno,'modi')
        all_tpm=f'{self.outdir}/all_merged.tpm'
        all_anno_tpm=general.insLable(all_tpm,'anno')
        all_anno_modi_tpm=general.insLable(all_tpm,'modi')
        cmd.extend(
            ["sed '1s/Sequence_ID/Contig/'",tax_anno,'>',modi_tax_anno,'\n',
            'linkTab.py',all_tpm,modi_tax_anno,'left Contig',all_anno_tpm,'\n',
            'tpmAddSource.py',all_anno_tpm,all_anno_modi_tpm,'\n',
            'pheatmap_for_tpm.R',all_anno_modi_tpm,self.samp_info,self.outdir,'\n']
        )
        len_sum_tpm_qual_xls=f'{self.outdir}/contig_quality_summary.xls'
        cmd.extend(
             ['sumAbundance.py',all_anno_modi_tpm,len_sum_tpm_qual_xls,self.outdir,'\n',
             'fa_length_tpm_scatter.R',len_sum_tpm_qual_xls,self.outdir,'\n']
        )
        tax_tpm=f'{self.outdir}/tax_tpm.xls'
        cmd.extend(
            ['abundByTax.py',all_anno_modi_tpm,self.outdir,'\n',
             'barplot_for_taxa_tpm.R',tax_tpm,self.outdir,'\n',
             'NMDS.R',all_anno_modi_tpm,self.samp_info,self.outdir,'\n']
        )
        shell=f'{self.outdir}/stat_plot.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
