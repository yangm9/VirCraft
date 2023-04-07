from .multiQuant import multiVirCount
from ..general import utils
from ..identify.viridsop import VirScan

class VirAbdStat(multiVirCount):
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def mergeAbd(self): #Heatmap for contigs abundance
        abd=f'{self.outdir}/all_merged.cov'
        cmd=['merge_abds.pl',self.samp_info,self.outdir,'cov\n']
        return cmd,abd
    def sizeAbdPlot(self,abd:str,checkv_dir=None):
        cmd=[]
        sum_abd=f'{self.outdir}/contig_sum_abd.xls'
        if not checkv_dir:
            tmp_fa=self.fasta
            tmp_dir=self.outdir
            tmp_threads=self.threads
            VirSOP=VirScan(tmp_fa,tmp_dir,tmp_threads)
            cmd,__=VirSOP.checkv(tmp_fa)
            checkv_dir=f'{self.outdir}/checkv'
        votu_qual=f'{checkv_dir}/quality_summary.tsv'
        sum_qual=f'{self.outdir}/sum_abd_qual.xls'
        cmd.extend(
            ['sum_abd_by_seq.py',abd,sum_abd,'\n',
            "sed -i '1s/contig_id/Contig/'",votu_qual,'\n',
            'linkTab.py',sum_abd,votu_qual,'left Contig',sum_qual,'\n',
            'variables_scatter.R',sum_qual,'contig_length~Total_Abundance~checkv_quality',self.outdir,'\n',
            'variables_scatter.R',sum_qual,'contig_length~contamination~checkv_quality',self.outdir,'\n',
            'variables_scatter.R',sum_qual,'contig_length~gene_count~checkv_quality',self.outdir,'\n']
        )
        return cmd
    def taxaAbd(self,abd:str,taxa_anno=None):
        if not taxa_anno:
            return ['pheatmap_for_abd.R',abd,self.samp_info,self.outdir,'\n']
        m_taxa_anno=f'{self.outdir}/DemoVir_assignments.txt'
        abd_taxa=f'{self.outdir}/all_sum_abd_taxa.xls'
        m_abd_taxa=f'{self.outdir}/all_sum_abd_taxa.m.xls'
        ctg_taxa_abd=f'{self.outdir}/all_ctg_abd_taxa.xls'
        taxa_sum_abd=f'{self.outdir}/all_taxa_sum_abd.xls'
        cmd=['echo "Heatmap for abundance"\n',
            "sed '1s/Sequence_ID/Contig/'",taxa_anno,'>',m_taxa_anno,'\n',
            'linkTab.py',abd,m_taxa_anno,'left Contig',abd_taxa,'\n',
            "sed '1s/Order/Source/'",abd_taxa,'>',m_abd_taxa,'\n',
            'pheatmap_for_abd.R',m_abd_taxa,self.samp_info,self.outdir,'\n',
            'echo "Barplot for abundance by taxa"\n',
            'taxa_annot_abd.py',abd_taxa,ctg_taxa_abd,'\n',
            'sum_abd_by_taxa.py',ctg_taxa_abd,self.outdir,'\n',
            'barplot_for_taxa_abd.R',taxa_sum_abd,self.outdir,'\n']
        return cmd
    def diversity(self,abd:str):
        alpha_diversity=f'{self.outdir}/alpha_diversity.xls'
        cmd=['echo "Alpha and Beta Diversity"\n',
            'alpha_diversity.R',abd,alpha_diversity,'\n',
            'NMDS.R',abd,self.samp_info,self.outdir,'\n']
        return cmd
    def QuantStat(self,taxa_anno=None,checkv_dir=None,batch_size=2):
        batch_size=str(batch_size)
        self.virCountBySamp()
        cmd=[self.envs]
        cmd.extend(
            ['multithreads.pl',self.outdir,'viral_count.sh',batch_size,'\n']
        )
        tmp_cmd,abd=self.mergeAbd()
        cmd.extend(tmp_cmd)
        cmd.extend(self.sizeAbdPlot(abd,checkv_dir))
        cmd.extend(self.taxaAbd(abd,taxa_anno))
        cmd.extend(self.diversity(abd))
        shell=f'{self.outdir}/{self.name}_vir_count.sh'
        utils.printSH(shell,cmd)
        #results=utils.execute(cmd)
        return results
