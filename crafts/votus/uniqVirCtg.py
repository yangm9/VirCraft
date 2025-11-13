from os import path
from ..general import utils
from ..identify.viralDetectors import VirDetectTools

class VirRef(VirDetectTools):
    '''
    Generate vOTUs
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir, threads)
    def cdhit_cluster(self):
        '''
        Cluster the sequence and remove redundancy for FastA file using CD-HIT.
        '''
        votus = f'{self.outdir}/{self.name}_votus.fa'
        cmd = [utils.selectENV('VC-General')]
        cmd = ['cd-hit-est', '-i', self.fasta, '-o', votus, '-T', self.threads, '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000\n']
        return cmd, votus
    def blast_cluster(self):
        '''
        Rapid genome clustering based on pairwise ANI provided by CheckV.
        '''
        blastdb = f'{self.wkfile_dir}/{self.name}.db'
        blast_out = f'{self.wkfile_dir}/{self.name}.blast'
        ani_out = f'{self.wkfile_dir}/{self.name}.ani'
        clusters = f'{self.wkfile_dir}/{self.name}.clusters'
        votu_list = f'{self.wkfile_dir}/{self.name}.votulist'
        votus = f'{self.wkfile_dir}/{self.name}_votus.fa'
        cmd = ['makeblastdb', '-in', self.fasta, '-dbtype nucl', '-out', blastdb, '\n',
               'blastn', '-query', self.fasta, '-db', blastdb, '-num_threads', self.threads, '-out', blast_out, "-outfmt '6 std qlen slen' -max_target_seqs 10000\n",
               'anicalc.py', '-i', blast_out, '-o', ani_out,'\n',
               'aniclust.py', '--fna', self.fasta, '--ani', ani_out, '--out', clusters, '--min_ani 95 --min_tcov 85 --min_qcov 0\n',
               'cut -f 1', clusters, '>', votu_list, '\n',
               'extrSeqByName.pl', votu_list, self.fasta, votus, '\n']
        cmd.extend(
           ['cd', self.outdir, '&& ln', ani_out, '&& ln', clusters, '&& ln', votus, '\n']
        )
        return cmd, votus
    def cluster(self, method):
        if method == 'blast':
            return self.blast_cluster()
        elif method == 'cdhit':
            return self.cdhit_cluster()
        else:
            raise Exception('method should only be "blast" or "cdhit".')
    def votuQC(self, votus, min_len=2000):
        cmd, __ = self.checkv(votus)
        checkv_qual = f'{self.wkfile_dir}/checkv/quality_summary.tsv'
        cmd.append(utils.selectENV('VC-General'))
        cmd.extend(
            ['pie_plot.R', checkv_qual, 'checkv_quality', self.stat_dir, '\n',
             'pie_plot.R', checkv_qual, 'provirus', self.stat_dir, '\n',
             'pie_plot.R', checkv_qual, 'miuvig_quality', self.stat_dir, '\n',
             'quality_length_boxplot.R', checkv_qual, self.stat_dir, '\n']
        )
        original_fasta = self.fasta
        self.fasta = votus
        cmd.extend(self.statFA())
        self.fasta = original_fasta #self.fasta need to be changed back to its original value
        original_outdir = self.outdir
        self.outdir = self.wkfile_dir
        tmp_cmd, vbdir = self.vibrant(str(min_len))
        self.outdir = original_outdir #self.outdir need to be changed back to its original value
        cmd.extend(tmp_cmd)
        votus_prefix = f'{self.name}_votus'
        vb_vir_info = f'{self.wkfile_dir}/VIBRANT_{votus_prefix}/VIBRANT_results_{votus_prefix}/VIBRANT_genome_quality_{votus_prefix}.tsv'
        vb_ckv_tsv = f'{self.stat_dir}/votus_lifetype_quality.tsv'
        cmd.append(utils.selectENV('VC-General'))
        cmd.extend(
            ['votus_lifetype_quality.py', checkv_qual, vb_vir_info, vb_ckv_tsv, '\n',
             'votus_lifetype_quality_barplot.R', vb_ckv_tsv, self.stat_dir, '\n']
        )
        return cmd
    def vCTGs_cluster(self, cov_input=None, checkv=None, unrun=False, method='blast', min_len=2000):
        cmd = []
        # Gene prediction using prodigal
        tmp_cmd, votus = self.cluster(method)
        cmd.extend(tmp_cmd)
        cmd.extend(self.votuQC(votus, min_len))
        shell = f'{self.shell_dir}/{self.name}_votu_cluster.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
