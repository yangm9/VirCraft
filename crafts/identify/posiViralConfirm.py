import os
from ..general import utils
from .viralDetectors import VirDetectTools

class vIdentify(VirDetectTools):
    '''
    Main Scripts of identify module.
    self.BATCH_SIZE initialized as 4 in VirCfg class will be used to divide the input threads into 2*self.BATCH_SIZE portions, with 2 allocated to VirSorter2, and the other 2 allocated to VIBRANT and DeepVirFinder respectively.
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.allthreads = threads
        self.threads = int(threads) // (self.BATCH_SIZE * 2)
    def vFilter(self, min_len=2000, filt_mode='permissive'):
        score_tsv = f'{self.wkfile_dir}/all_viral_ctgs.score.tsv'
        score_filt_tsv = utils.insLable(score_tsv, filt_mode)
        viral_filt_ctg_list = f'{self.wkfile_dir}/viral_filt_ctg.list'
        viral_filt_ctgs_fna = f'{self.wkfile_dir}/viral_filt_ctg.fna'
        viral_posi_ctgs_fna = f'{self.wkfile_dir}/viral_positive_ctg.fna'
        vs2_list = f'{self.wkfile_dir}/VirSorter2.contigs'
        vb_list = f'{self.wkfile_dir}/VIBRANT.contigs'
        dvf_list = f'{self.wkfile_dir}/DeepVirFinder.contigs'
        gn_list = f'{self.wkfile_dir}/geNomad.contigs'
        venn_pdf = f'{self.stat_dir}/viral_contigs_multitools_comparison_venn.pdf'
        cmd = [utils.selectENV('VC-General')]
        tmp_cmd = ''
        #tmp_cmd, cat_dir = self.contig_annotation_tool(viral_filt_ctgs_fna)
        cmd.extend(
            ['merge_vctg_info.py', self.name, self.wkfile_dir, filt_mode, # generate all_viral_ctgs.score.tsv file 
             '&& cut -f 1', score_filt_tsv, "|sed '1d' >", viral_filt_ctg_list, '&& extrSeqByName.pl', viral_filt_ctg_list, self.fasta, viral_filt_ctgs_fna, '\n\n',
             "awk -F '\\t' 'NR>1 && $32>=1 {print $1}'", score_tsv, '>', vs2_list, "&& awk -F '\\t' 'NR>1 && $33>=1 {print $1}'", score_tsv, '>', vb_list,
             "&& awk -F '\\t' 'NR>1 && $34>=1 {print $1}'", score_tsv, '>', dvf_list, "&& awk -F '\\t' 'NR>1 && $35>=1 {print $1}'", score_tsv, '>', gn_list,
             '&& venn4.R', vs2_list, vb_list, dvf_list, gn_list, venn_pdf, '\n\n']
        )
        tmp_cmd, checkv_dir = self.checkv(viral_filt_ctgs_fna)
        cmd.extend(tmp_cmd)
        quality_summary_tsv = checkv_dir + '/quality_summary.tsv'
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(
            ['vir_qual_filt.py', score_filt_tsv, checkv_dir, viral_filt_ctgs_fna, viral_posi_ctgs_fna, '\n']
        )
        original_fasta = self.fasta
        self.fasta = viral_posi_ctgs_fna
        cmd.extend(self.statFA()) # Invoke the statFA() method from the Seq module
        self.fasta = original_fasta
        complete_ctg_fna = viral_posi_ctgs_fna.replace('.fna', '_complete.fna')
        provirus_ctg_fna = viral_posi_ctgs_fna.replace('.fna', '_provirus.fna')
        vctg_for_binning_fna = viral_posi_ctgs_fna.replace('.fna', '_for_binning.fna')
        cmd.extend(
            ['cp', score_tsv, self.outdir, '&& cp', viral_posi_ctgs_fna, self.outdir,
             '&& cp', complete_ctg_fna, self.outdir, '&& cp', provirus_ctg_fna, self.outdir,
             '&& cp', vctg_for_binning_fna, self.outdir, '\n']
        )
        return cmd
    def Identify(self, min_len=2000, filt_mode='permissive', unrun=False):
        try:
            if int(self.allthreads) < 8:
                raise ValueError('The threads number must not be less than 8!!!')
        except ValueError as e:
            print(f'ERROR: {e}')
            exit(1)
        self.threads = str(self.threads)
        #vibrant
        cmd, wkdir = self.vibrant(min_len)
        shell = f'{self.shell_dir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell, cmd)
        #deepvirfinder
        cmd, wkdir = self.deepvirfinder(min_len)
        shell = f'{self.shell_dir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell, cmd)
        #genomad
        cmd, wkdir = self.genomad()
        shell = f'{self.shell_dir}/{self.name}_gn_ctg.sh'
        utils.printSH(shell, cmd)
        #VirSorter2
        self.threads = str(int(self.allthreads) - (int(self.threads) * 3)) #5
        cmd, wkdir = self.virsorter(in_fa=self.fasta, n=0, min_len=min_len, min_score=0.5)
        shell = f'{self.shell_dir}/{self.name}_vs2_ctg.sh'
        utils.printSH(shell, cmd)
        #multiple run
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['multithreads.pl', self.shell_dir, 'ctg.sh 4\n']
        )
        shell = f'{self.shell_dir}/{self.name}_batch_identify_virus.sh'
        utils.printSH(shell, cmd)
        results = ''
        if not unrun: results = utils.execute(shell) 
        self.threads = str(self.allthreads) #8
        cmd = self.vFilter(min_len, filt_mode)
        shell = f'{self.shell_dir}/{self.name}_get_positive_virus.sh'
        utils.printSH(shell, cmd)
        if not unrun: results += utils.execute(shell)
        return results
