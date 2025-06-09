import os
from ..general import utils
from ..data.bioseq import Seq
from ..identify.viralDetectors import VirDetectTools

class vIdentify(VirDetectTools):
    '''
    Main Scripts of identify module.
    self.BATCH_SIZE initialized as 4 in VirCfg class will be used to divide the input threads into 2*self.BATCH_SIZE portions, with 2 allocated to VirSorter2, and the other 2 allocated to VIBRANT and DeepVirFinder respectively.
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.allthreads = threads
        self.threads = int(threads) // (self.BATCH_SIZE * 2)
    def vFilter(self, cutoff=1500, mode='permissive'):
        mode_dict = {'permissive' : '1', 'strict' : '2'}
        score_tsv = f'{self.wkfile_dir}/all_viral_ctgs.score.tsv'
        score_filt_tsv = utils.insLable(score_tsv, mode)
        viral_filt_ctg_list = f'{self.wkfile_dir}/viral_filt_ctg.list'
        viral_filt_ctgs_fna = f'{self.wkfile_dir}/viral_filt_ctg.fna'
        viral_posi_ctgs_fna = f'{self.wkfile_dir}/viral_positive_ctg.fna'
        vs2_list = f'{self.wkfile_dir}/vs2_ctg.list'
        vb_list = f'{self.wkfile_dir}/vb_ctg.list'
        dvf_list = f'{self.wkfile_dir}/dvf_ctg.list'
        gn_list = f'{self.wkfile_dir}/gn_ctg.list'
        venn_pdf = f'{self.wkfile_dir}/vs2_vb_dvf_gn_venn.pdf'
        cmd = [utils.selectENV('VC-General')]
        tmp_cmd = ''
        #tmp_cmd, cat_dir = self.contig_annotation_tool(viral_filt_ctgs_fna)
        cmd.extend(
            ['merge_ctg_list.py', self.name, self.wkfile_dir, f"&& awk -F '\\t' 'NR == 1 || $32 >= {mode_dict[mode]}'", score_tsv, '>', score_filt_tsv, '\n',
            'cut -f 1', score_filt_tsv, "|sed '1d' >", viral_filt_ctg_list, '&& extrSeqByName.pl', viral_filt_ctg_list, self.fasta, viral_filt_ctgs_fna, '\n',
            "awk -F '\\t' '$28>=1 {print $1}'", score_tsv, '>', vs2_list, "&& awk -F '\\t' '$29>=1 {print $1}'", score_tsv, '>', vb_list,
            "&& awk -F '\\t' '$30>=1 {print $1}'", score_tsv, '>', dvf_list, "&& awk -F '\\t' '$31>=1 {print $1}'", score_tsv, '>', gn_list,
            '&& venn4.R', vs2_list, vb_list, dvf_list, gn_list, venn_pdf, '\n']
        )
        if cutoff <= 5000:
            tmp_cmd, checkv_dir = self.checkv(viral_filt_ctgs_fna)
            cmd.extend(tmp_cmd)
            quality_summary_tsv = checkv_dir + '/quality_summary.tsv'
            cmd.extend([utils.selectENV('VC-General')])
            cmd.extend(
                ['vir_qual_filt.py', quality_summary_tsv, viral_filt_ctgs_fna, viral_posi_ctgs_fna, '\n']
            )
        else:
            cmd.extend(
                ['cp', viral_filt_ctgs_fna, viral_posi_ctgs_fna, '\n']
            )
        outdir = self.outdir
        self.outdir = viral_posi_ctgs_fna
        cmd.extend(self.statFA()) # Invoke the statFA() method from the Seq module
        self.outdir = outdir
        cmd.extend(
            ['cd', self.outdir, '&& ln', score_tsv, '&& ln', viral_posi_ctgs_fna, '\n']
        )
        return cmd
    def Identify(self, cutoff=1500, mode='permissive', unrun=False):
        try:
            if int(self.allthreads) < 8:
                raise ValueError('The threads number must not be less than 8!!!')
        except ValueError as e:
            print(f'ERROR: {e}')
            exit(1)
        self.threads = str(self.threads)
        #vibrant
        cmd, wkdir = self.vibrant(cutoff)
        shell = f'{self.shell_dir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell, cmd)
        #deepvirfinder
        cmd, wkdir = self.deepvirfinder(cutoff)
        shell = f'{self.shell_dir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell, cmd)
        #genomad
        cmd, wkdir = self.genomad()
        shell = f'{self.shell_dir}/{self.name}_gm_ctg.sh'
        utils.printSH(shell, cmd)
        #VirSorter2
        self.threads = str(int(self.allthreads) - (int(self.threads) * 3)) #5
        cmd, wkdir = self.virsorter(in_fa=self.fasta, n=0, min_length=cutoff, min_score=0.5)
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
        cmd = self.vFilter(cutoff, mode)
        shell = f'{self.shell_dir}/{self.name}_get_positive_virus.sh'
        utils.printSH(shell, cmd)
        if not unrun: results += utils.execute(shell)
        return results
