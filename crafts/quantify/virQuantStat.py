import os
from .multiQuant import multiVirCount
from ..general import utils
from ..identify.viralDetectors import VirDetectTools

#Main Class
class VirAbdStat(multiVirCount):
    def __init__(self, samp_info=None, fasta=None, outdir=None, threads=8):
        super().__init__(samp_info, fasta, outdir, threads)
    def mergeAbd(self, seq_type='Contig'): #Heatmap for contigs abundance
        abd = f'{self.stat_dir}/all_viral_contig_abundance.cov'
        cmd = ['echo "Merge all viral contig abundance files"\n'
               'paste_abundance_files.py', self.wkfile_dir, self.samp_info, abd, seq_type, '\n\n'
        ]
        return cmd, abd
    def sizeAbdPlot(self, abd: str, checkv_dir=None):
        cmd = ['echo "Draw scatter plots of abundance and quality information"\n']
        if not checkv_dir:
            tmp_fa = self.fasta
            tmp_threads = self.threads
            VirDetect = VirDetectTools(fasta=tmp_fa, outdir=self.outdir, threads=tmp_threads)
            tmp_cmd, __ = VirDetect.checkv(tmp_fa)
            checkv_dir = f'{self.wkfile_dir}/checkv'
            cmd.extend(tmp_cmd)
        checkv_dir = os.path.abspath(checkv_dir)
        cmd.extend([utils.selectENV('VC-General')])
        ctg_sum_abd_tsv = f'{self.stat_dir}/viral_contig_total_abundance.tsv'
        votu_qual_tsv = f'{checkv_dir}/quality_summary.tsv'
        ctg_sum_qual_tsv = f'{self.stat_dir}/vctg_total_abundance_quality.tsv'
        cmd.extend(
            ['sum_abd_by_seq.py', abd, ctg_sum_abd_tsv, "&& sed -i '1s/contig_id/Contig/'", votu_qual_tsv,
             '&& linkTab.py', ctg_sum_abd_tsv, votu_qual_tsv, 'left Contig', ctg_sum_qual_tsv,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~Total_Abundance~checkv_quality', self.stat_dir,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~contamination~checkv_quality', self.stat_dir,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~gene_count~checkv_quality', self.stat_dir, '\n\n']
        )
        return cmd
    def taxaAbd(self, abd: str, taxa_anno=None):
        cmd = [utils.selectENV('VC-General')]
        if not taxa_anno:
            cmd.extend(
                ['pheatmap_for_abd.R', abd, self.samp_info, self.stat_dir, '\n']
            )
            return cmd
        taxa_anno = os.path.abspath(taxa_anno)
        vir_taxa_anno_tsv = f'{self.stat_dir}/vctg_taxa.tsv'
        vir_abd_taxa_tsv = f'{self.stat_dir}/vctg_abundance_order_family.tsv'
        vir_taxa_abd_tsv = f'{self.stat_dir}/vctg_taxa_abundance.tsv'
        vir_abd_sum_by_taxa_tsv = f'{self.stat_dir}/summed_vtaxa_abundance.tsv'
        order_sum_vir_abd_tsv = f'{self.stat_dir}/summed_vOrder_abundance.tsv'
        family_sum_vir_abd_tsv = f'{self.stat_dir}/summed_vFamily_abundance.tsv'
        cmd.extend(
            ['echo "Ploting heatmaps of relative abundance"\n',
             "sed '1s/Sequence_ID/Contig/'", taxa_anno, '>', vir_taxa_anno_tsv,
             '&& linkTab.py', abd, vir_taxa_anno_tsv, 'left Contig', vir_abd_taxa_tsv,
             '&& pheatmap_for_abd.R', vir_abd_taxa_tsv, self.samp_info, self.stat_dir, 'Order',
             '&& pheatmap_for_abd.R', vir_abd_taxa_tsv, self.samp_info, self.stat_dir, 'Family\n\n',
             
             'echo "Barplot for abundance by taxa"\n',
             'taxa_annot_abd.py', vir_abd_taxa_tsv, vir_taxa_abd_tsv,
             '&& sum_abd_by_taxa.py', vir_taxa_abd_tsv, self.stat_dir, '\n',
             'barplot_for_taxa_abd.R', vir_abd_sum_by_taxa_tsv, self.stat_dir,
             '&& barplot_for_taxa_abd.R', order_sum_vir_abd_tsv, self.stat_dir,
             '&& barplot_for_taxa_abd.R', family_sum_vir_abd_tsv, self.stat_dir, '\n\n',
             
             'echo "Heatmap for abundance by taxa"\n',
             'pheatmap_for_abd.R', order_sum_vir_abd_tsv, self.samp_info, self.stat_dir,
             '&& pheatmap_for_abd.R', family_sum_vir_abd_tsv, self.samp_info, self.stat_dir, '\n\n']
        )
        return cmd
    def diversity(self, abd: str):
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['echo "Alpha and Beta diversity analysis"\n',
             'alpha_diversity.R', abd, self.samp_info, self.stat_dir, '\n',
             'PCoA.R', abd, self.samp_info, self.stat_dir, '\n',
             'NMDS.R', abd, self.samp_info, self.stat_dir, '\n\n']
        )
        return cmd
    def QuantStat(self, taxa_anno=None, checkv_dir=None, coverm_method='mean', unrun=False, clear=False): # Main Function
        cmd = self.virCountBySamp(coverm_method)
        cmd.append('echo "Run coverm in batch mode to calculate the abundance"\n')
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(
            ['multithreads.pl', self.shell_dir, 'viral_count.sh', str(self.BATCH_SIZE), '\n\n']
        )
        tmp_cmd, abd = (None, None)
        if coverm_method == 'metabat':
            tmp_cmd, abd = self.mergeAbd(coverm_method) 
            cmd.extend(tmp_cmd)
        else:
            tmp_cmd, abd = self.mergeAbd()
            cmd.extend(tmp_cmd)
            cmd.extend(self.sizeAbdPlot(abd, checkv_dir))
            cmd.extend(self.taxaAbd(abd, taxa_anno))
            cmd.extend(self.diversity(abd))
        shell = f'{self.shell_dir}/{self.name}_vir_count.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
