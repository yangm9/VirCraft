from ..general import cmdExec,general
from .align import multiAln

class VirTPM(multiAln):
    '''
    '''
    def __init__(self,fastqs='',outdir='',threads=8):
        super().__init__(self,config,outdir,threads)
    def calcTPM(self,samp:str,sort_bam:str,wkdir:str):
        cmd=[self.envs]
        tpm=f'{wkdir}/{samp}.tpm'
        cmd.extend(
            ['coverm contig','-b',sort_bam,
             '-t 20 --min-read-aligned-length 50 --min-read-percent-identity 0.95 --proper-pairs-only -m tpm',
             '>',tpm]
        )
        shell=f'{wkdir}/{samp}_tpm.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def tpmBySamp(self):
        results=''
        wkdir=f'{self.wkdir}/2.coverm'
        general.mkdir(wkdir)
        bam_dir=f'{self.wkdir}/1.bwa'
        for samp in self.sampDict.keys():
            sort_bam=f'{bam_dir}/{samp}.sort.bam'
            results+=self.calcTPM(samp,sort_bam,wkdir)
        cmd=['merge_tpms.pl',self.confDict['SampInfo'],wkdir,'\n']
        modi_tax_anno=general.insLable(self.tax_anno,'modi')
        merged_tpm=f'{wkdir}/all_merged.tpm'
        merged_anno_tpm=general.insLable(merged_tpm,'anno')
        merged_anno_modi_tpm=general.insLable(merged_tpm,'modi')
        cmd.extend(
            ["sed '1s/Sequence_ID/Contig/'",self.tax_anno,'>',modi_tax_anno,'\n',
             'linkTab.py',merged_tpm,modi_tax_anno,'left Contig',merged_anno_tpm,'\n',
             'tpmAddSource.py',merged_anno_tpm,merged_anno_modi_tpm,'\n',
             'pheatmap_for_tpm.R',merged_anno_modi_tpm,confDict['SampInfo'],wkdir,'\n']
        )
        len_sum_tpm_qual_xls=f'{wkdir}/contig_quality_summary.xls'
        cmd.extend(
             ['sumAbundance.py',merged_anno_modi_tpm,self.qual_summ,wkdir,'\n','fa_length_tpm_scatter.R',len_sum_tpm_qual_xls,wkdir,'\n']
        )
        tax_tpm=f'{wkdir}/tax_tpm.xls'
        cmd.extend(
            ['abundByTax.py',merged_anno_modi_tpm,wkdir,'\n',
             'barplot_for_taxa_tpm.R',tax_tpm,wkdir,'\n']
        )
        shell=f'{wkdir}/merged_anno_tpm.sh'
        general.printSH(shell,cmd)
        results+=cmdExec.execute(cmd)
        return results
