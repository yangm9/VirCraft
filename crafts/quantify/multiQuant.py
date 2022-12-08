import os
import glob
from ..general import cmdExec,general
from .alnQuant import VirCount
from ..config.config import Seq

class multiCount(Seq):
    def __init__(self,fastqs='',fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.fastqs=glob.glob(os.path.abspath(fastqs))
        self.indir=os.path.dirname(self.fastqs[0])
        self.postfix=os.path.splitext(self.fastqs[0])[1]
        fq1postfix=f'_1{self.postfix}'
        self.samps=[os.path.basename(name.replace(fq1postfix,'')) for name in fastqs if name.endswith(fq1postfix)]
        self.samp_list_f=f'{self.outdir}/samp.list'
        self.threads=str(threads)
    def countBySamp(self,group='NA'):
        _,bwa_idx=self.mkBwaIdx
        results=''
        wkdir=f'{self.outdir}/quant'
        general.mkdir(wkdir)
        SampList=open(self.samp_list_f,'w')
        for samp in self.samps:
            cmd=[self.envs]
            fq1=f'{self.indir}/{samp}_1{self.postfix}'
            fq2=f'{self.indir}/{samp}_2{self.postfix}'
            Count=VirCount(fq1,fq2,wkdir,self.threads)
            cmd.extend(Count.aln(samp,bwa_idx))
            cmd.extend(Count.coverm(samp))
            SampList.write(f'{samp}\t{group}i\n')
            shell=f'{wkdir}/{samp}_count.sh'
            general.printSH(shell,cmd)
            results+=cmdExec.execute(cmd)
        return results
    def statPlot(self,tax_anno=None):
        cmd=[self.envs]
        wkdir=f'{self.outdir}/statistics'
        general.mkdir(wkdir)
        cmd.extend(['merge_tpms.pl',self.samp_list_f,wkdir,'\n'])
        modi_tax_anno=general.insLable(tax_anno,'modi')
        all_tpm=f'{wkdir}/all_merged.tpm'
        all_anno_tpm=general.insLable(merged_tpm,'anno')
        all_anno_modi_tpm=general.insLable(merged_tpm,'modi')
        cmd.extend(
            ["sed '1s/Sequence_ID/Contig/'",tax_anno,'>',modi_tax_anno,'\n',
            'linkTab.py',all_tpm,modi_tax_anno,'left Contig',all_anno_tpm,'\n',
            'tpmAddSource.py',all_anno_tpm,all_anno_modi_tpm,'\n',
            'pheatmap_for_tpm.R',all_anno_modi_tpm,self.samp_list_f,wkdir,'\n']
        )
        len_sum_tpm_qual_xls=f'{wkdir}/contig_quality_summary.xls'
        cmd.extend(
             ['sumAbundance.py',merged_anno_modi_tpm,self.qual_summ,wkdir,'\n',
             'fa_length_tpm_scatter.R',len_sum_tpm_qual_xls,wkdir,'\n']
        )
        tax_tpm=f'{wkdir}/tax_tpm.xls'
        cmd.extend(
            ['abundByTax.py',merged_anno_modi_tpm,wkdir,'\n',
             'barplot_for_taxa_tpm.R',tax_tpm,wkdir,'\n']
        )
        shell=f'{wkdir}/stat_plot.sh'
        general.printSH(shell,cmd)
        results+=cmdExec.execute(cmd)
        return results
