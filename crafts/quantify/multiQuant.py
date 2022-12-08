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
        _,bwa_idx=self.mkBwaIdx
        results=''
        wkdir=f'{self.outdir}/quant'
        general.mkdir(wkdir)
        for samp in self.sampDict.keys():
            cmd=[self.envs]
            fq1,fq2=self.sampDict[samp][1].split(',')
            Count=VirCount(fq1,fq2,wkdir,self.threads)
            cmd.extend(Count.aln(samp,bwa_idx))
            cmd.extend(Count.coverm(samp))
            shell=f'{wkdir}/{samp}_count.sh'
            general.printSH(shell,cmd)
            results+=cmdExec.execute(cmd)
        return results
    def statPlot(self,tax_anno='test.xls',qual_summ='test.txt'):
        cmd=[self.envs]
        wkdir=f'{self.outdir}/statistics'
        general.mkdir(wkdir)
        cmd.extend(['merge_tpms.pl',self.samp_info,wkdir,'\n'])
        modi_tax_anno=general.insLable(tax_anno,'modi')
        all_tpm=f'{wkdir}/all_merged.tpm'
        all_anno_tpm=general.insLable(all_tpm,'anno')
        all_anno_modi_tpm=general.insLable(all_tpm,'modi')
        cmd.extend(
            ["sed '1s/Sequence_ID/Contig/'",tax_anno,'>',modi_tax_anno,'\n',
            'linkTab.py',all_tpm,modi_tax_anno,'left Contig',all_anno_tpm,'\n',
            'tpmAddSource.py',all_anno_tpm,all_anno_modi_tpm,'\n',
            'pheatmap_for_tpm.R',all_anno_modi_tpm,self.samp_info,wkdir,'\n']
        )
        len_sum_tpm_qual_xls=f'{wkdir}/contig_quality_summary.xls'
        cmd.extend(
             ['sumAbundance.py',all_anno_modi_tpm,qual_summ,wkdir,'\n',
             'fa_length_tpm_scatter.R',len_sum_tpm_qual_xls,wkdir,'\n']
        )
        tax_tpm=f'{wkdir}/tax_tpm.xls'
        cmd.extend(
            ['abundByTax.py',all_anno_modi_tpm,wkdir,'\n',
             'barplot_for_taxa_tpm.R',tax_tpm,wkdir,'\n']
        )
        shell=f'{wkdir}/stat_plot.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
