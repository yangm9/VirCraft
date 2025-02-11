import os
from .multiQuant import multiVirCount
from ..general import utils
from ..identify.viridsop import VirScan

#Main Class
class VirAbdStat(multiVirCount):
    def __init__(self, samp_info=None, fasta=None, outdir=None,threads=8):
        super().__init__(samp_info, fasta, outdir, threads)
    def mergeAbd(self): #Heatmap for contigs abundance
        abd = f'{self.wkdir}/all_merged.cov'
        cmd = ['paste_abundance_files.py', self.samp_info, self.wkdir, 'cov\n']
        return cmd, abd
    def sizeAbdPlot(self, abd: str, checkv_dir=None):
        cmd = [utils.selectENV('VC-General')]
        sum_abd = f'{self.wkdir}/contig_sum_abd.tsv'
        if not checkv_dir:
            tmp_fa = self.fasta
            tmp_dir = self.wkdir
            tmp_threads = self.threads
            VirSOP = VirScan(tmp_fa, tmp_dir, tmp_threads)
            cmd,__ = VirSOP.checkv(tmp_fa)
            checkv_dir = f'{self.wkdir}/checkv'
        checkv_dir = os.path.abspath(checkv_dir)
        votu_qual = f'{checkv_dir}/quality_summary.tsv'
        sum_qual = f'{self.wkdir}/sum_abd_qual.tsv'
        cmd.extend(
            ['sum_abd_by_seq.py', abd,sum_abd, '\n',
            "sed -i '1s/contig_id/Contig/'", votu_qual, '\n',
            'linkTab.py', sum_abd, votu_qual, 'left Contig', sum_qual, '\n',
            'variables_scatter.R', sum_qual, 'contig_length~Total_Abundance~checkv_quality', self.wkdir, '\n',
            'variables_scatter.R', sum_qual, 'contig_length~contamination~checkv_quality', self.wkdir, '\n',
            'variables_scatter.R', sum_qual, 'contig_length~gene_count~checkv_quality', self.wkdir, '\n']
        )
        return cmd
    def taxaAbd(self, abd: str, taxa_anno=None):
        cmd = [utils.selectENV('VC-General')]
        if not taxa_anno:
            cmd.extend(
                ['pheatmap_for_abd.R', abd, self.samp_info, self.wkdir, '\n']
            )
            return cmd
        taxa_anno = os.path.abspath(taxa_anno)
        m_taxa_anno = f'{self.wkdir}/all_votu.taxa.txt'
        abd_taxa = f'{self.wkdir}/all_sum_abd_taxa.tsv'
        m_abd_taxa = f'{self.wkdir}/all_sum_abd_taxa.m.tsv'
        ctg_taxa_abd = f'{self.wkdir}/all_ctg_abd_taxa.tsv'
        taxa_sum_abd = f'{self.wkdir}/all_taxa_sum_abd.tsv'
        order_sum_abd = f'{self.wkdir}/all_Order_sum_abd.tsv'
        family_sum_abd = f'{self.wkdir}/all_Family_sum_abd.tsv'
        cmd.extend(
            ['echo "Heatmap for abundance"\n',
            "sed '1s/Sequence_ID/Contig/'", taxa_anno, '>', m_taxa_anno, '\n',
            'linkTab.py', abd, m_taxa_anno, 'left Contig', abd_taxa, '\n',
            "sed '1s/Order/Source/'", abd_taxa, '>', m_abd_taxa, '\n',
            'pheatmap_for_abd.R', m_abd_taxa, self.samp_info, self.wkdir, '\n',
            'echo "Barplot for abundance by taxa"\n',
            'taxa_annot_abd.py', abd_taxa, ctg_taxa_abd, '\n',
            'sum_abd_by_taxa.py', ctg_taxa_abd, self.wkdir, '\n',
            'barplot_for_taxa_abd.R', taxa_sum_abd, self.wkdir, '\n',
            'barplot_for_taxa_abd.R', order_sum_abd, self.wkdir, '\n',
            'barplot_for_taxa_abd.R', family_sum_abd, self.wkdir, '\n']
        )
        return cmd
    def diversity(self, abd: str):
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['echo "Alpha and Beta Diversity"\n',
            'alpha_diversity.R', abd, self.samp_info, self.wkdir, '\n',
            'NMDS.R', abd, self.samp_info, self.wkdir,'\n']
        )
        return cmd
    def QuantStat(self, taxa_anno=None, checkv_dir=None, unrun=False, clear=False):
        "Main Function"
        cmd = self.virCountBySamp()
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(
            ['multithreads.pl', self.shelldir, 'viral_count.sh', str(self.BATCH_SIZE), '\n']
        )
        tmp_cmd, abd = self.mergeAbd()
        cmd.extend(tmp_cmd)
        cmd.extend(self.sizeAbdPlot(abd, checkv_dir))
        cmd.extend(self.taxaAbd(abd, taxa_anno))
        cmd.extend(self.diversity(abd))
        shell = f'{self.shelldir}/{self.name}_vir_count.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
