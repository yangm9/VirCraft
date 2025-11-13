import os
from .multiQuant import multiVirCount
from ..general import utils
from ..identify.viralDetectors import VirDetectTools

# VirAbdStat: main class for viral diversity analysis based on viral contig abundances.
# This class is a tool for statistical analysis of viral abundance, with main functions including: 
# 1) Abundance calculation: Using the coverm tool to calculate the abundance of viral contigs in each sample.
# 2) Data integration: Merging the abundance data of all samples. Quality assessment: Analyzing the relationship between contig quality and abundance in combination with CheckV results. Diversity analysis: Alpha and Beta diversity analysis (PCoA, NMDS). 
# 3) Taxonomic analysis: Summarizing abundance by taxonomic levels (order, family, etc.) and visualizing. 
# 4) Visualization: Generating various charts such as scatter plots, bar charts, and heat maps.
class VirAbdStat(multiVirCount):
    def __init__(self, samp_info=None, fasta=None, outdir=None, threads=8):
        super().__init__(samp_info, fasta, outdir, threads)
        self.all_vctg_abd_cov = f'{self.stat_dir}/all_viral_contig_abundance.cov'
    def mergeAbd(self, seq_type='Contig'): #Heatmap for contigs abundance
        cmd = ['# Merge all viral contig abundance files\n'
               'paste_abundance_files.py', self.wkfile_dir, self.samp_info, self.all_vctg_abd_cov, seq_type, '\n\n'
        ]
        return cmd
    def sizeAbdPlot(self, checkv_dir=None):
        cmd = ['# Draw scatter plots of abundance and quality information\n']
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
            ['sum_abd_by_seq.py', self.all_vctg_abd_cov, ctg_sum_abd_tsv, "&& sed -i '1s/contig_id/Contig/'", votu_qual_tsv,
             '&& linkTab.py', ctg_sum_abd_tsv, votu_qual_tsv, 'left Contig', ctg_sum_qual_tsv,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~Total_Abundance~checkv_quality', self.stat_dir,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~contamination~checkv_quality', self.stat_dir,
             '&& variables_scatter.R', ctg_sum_qual_tsv, 'contig_length~gene_count~checkv_quality', self.stat_dir, '\n\n']
        )
        return cmd
    def diversity(self):
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['# Alpha and Beta diversity analysis"\n',
             'alpha_diversity.R', self.all_vctg_abd_cov, self.samp_info, self.stat_dir, '\n',
             'PCoA.R', self.all_vctg_abd_cov, self.samp_info, self.stat_dir, '\n',
             'NMDS.R', self.all_vctg_abd_cov, self.samp_info, self.stat_dir, '\n\n']
        )
        return cmd
    def taxaAbdPlot(self, taxa_anno=None):
        cmd = [utils.selectENV('VC-General')]
        #if not taxa_anno:
        #    cmd.extend(
        #        ['pheatmap_for_abd.R', self.all_vctg_abd_cov, self.samp_info, self.stat_dir, '\n']
        #    )
        #    return cmd
        taxa_anno = os.path.abspath(taxa_anno)
        vir_taxa_anno_tsv = f'{self.stat_dir}/vctg_taxa.tsv' 
        vir_abd_taxa_tsv = f'{self.stat_dir}/vctg_abundance_taxa.tsv' # Columns: Contig,Sample_1,Sample_2...Sample_N,Superrealm,...,Species stored abundance information
        vir_abd_sum_by_taxa_tsv = f'{self.stat_dir}/summed_vtaxa_abundance.tsv'
        order_sum_vir_abd_tsv = f'{self.stat_dir}/summed_vOrder_abundance.tsv'
        family_sum_vir_abd_tsv = f'{self.stat_dir}/summed_vFamily_abundance.tsv'
        cmd.extend(
            ['# Barplot for abundance by taxa\n',
             "sed '1s/Sequence_ID/Contig/'", taxa_anno, '>', vir_taxa_anno_tsv,
             '&& linkTab.py', self.all_vctg_abd_cov, vir_taxa_anno_tsv, 'left Contig', vir_abd_taxa_tsv,
             #'&& pheatmap_for_abd.R', vir_abd_taxa_tsv, self.samp_info, self.stat_dir, 'Order',
             #'&& pheatmap_for_abd.R', vir_abd_taxa_tsv, self.samp_info, self.stat_dir, 'Family\n\n',
             '&& sum_abd_by_taxa.py', vir_abd_taxa_tsv, self.stat_dir, '\n',
             'barplot_for_taxa_abd.R', vir_abd_sum_by_taxa_tsv, self.stat_dir,
             '&& barplot_for_taxa_abd.R', order_sum_vir_abd_tsv, self.stat_dir,
             '&& barplot_for_taxa_abd.R', family_sum_vir_abd_tsv, self.stat_dir, '\n\n',
             '# Heatmap for abundance by taxa\n',
             'pheatmap_for_abd.R', order_sum_vir_abd_tsv, self.samp_info, self.stat_dir,
             '&& pheatmap_for_abd.R', family_sum_vir_abd_tsv, self.samp_info, self.stat_dir, '\n\n']
        )
        return cmd
    def QuantStat(self, votu_table=None, taxa_anno=None, checkv_dir=None, coverm_method='mean', unrun=False, clear=False): # Main Function
        cmd = ['# Run BWA, samtools and CoverM in batch mode to calculate the abundance\n']
        if votu_table:
            cmd.extend(['# Using provided votu_table, skipping abundance calculation steps\n'])
            cmd.extend(['cp', votu_table, self.all_vctg_abd_cov, '\n'])
        else:
            cmd.extend(self.virCountBySamp(coverm_method))
            cmd.extend([utils.selectENV('VC-General')])

            cmd.extend(
                ['multithreads.pl', self.shell_dir, 'viral_count.sh', str(self.BATCH_SIZE), '\n\n']
            )
            tmp_cmd = None
            if coverm_method == 'metabat':
                tmp_cmd = self.mergeAbd(coverm_method) 
            else:
                tmp_cmd = self.mergeAbd()
            cmd.extend(tmp_cmd)
        if not votu_table or coverm_method != 'metabat': # The subsequent steps are executed regardless of whether there is a votu_table or not.
            cmd.extend(self.sizeAbdPlot(checkv_dir))
            cmd.extend(self.diversity())
            cmd.extend(self.taxaAbdPlot(taxa_anno))
        shell = f'{self.shell_dir}/{self.name}_vir_count.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
