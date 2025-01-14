import os
import sys
from ..general import utils
from ..identify.multiFind import MultiTools

class vIdentify(MultiTools):
    '''
    Main Scripts of identify module.
    self.BATCH_SIZE initialized as 4 in VirCfg class will be used to divide the input threads into 2*self.BATCH_SIZE portions, with 2 allocated to VirSorter2, and the other 2 allocated to VIBRANT and DeepVirFinder respectively.
    '''
    def __init__(self,fasta = '', outdir = '', threads = 8):
        super().__init__(fasta, outdir)
        self.allthreads = threads
        self.threads = int(threads) // (self.BATCH_SIZE * 2)
    def vFilter(self, cutoff = 1500):
        score_tsv = f'{self.outdir}/all_viral_ctgs.score.tsv'
        score_filt_tsv = utils.insLable(score_tsv, 'gt2')
        viral_filt_ctg_list = f'{self.outdir}/viral_filt_ctg.list'
        viral_filt_ctgs_fna = f'{self.outdir}/viral_filt_ctg.fna'
        viral_posi_ctgs_fna = f'{self.outdir}/viral_positive_ctg.fna'
        cmd = [utils.selectENV('VC-General')]
        tmp_cmd = ''
        #tmp_cmd, cat_dir = self.contig_annotation_tool(viral_filt_ctgs_fna)
        cmd.extend(
            ['merge_ctg_list.py', self.name, self.outdir, "&& awk -F '\\t' 'NR==1 || $21>=2'", score_tsv, '>', score_filt_tsv, '\n',
            'cut -f 1', score_filt_tsv, "|sed '1d' >", viral_filt_ctg_list, '&& extrSeqByName.pl', viral_filt_ctg_list, self.fasta, viral_filt_ctgs_fna, '\n']
        )
        if cutoff <= 5000:
            tmp_cmd, checkv_fa = self.checkv(viral_filt_ctgs_fna)
            cmd.extend(tmp_cmd)
            quality_summary_tsv = os.path.dirname(checkv_fa) + '/quality_summary.tsv'
            cmd.extend([utils.selectENV('VC-General')])
            cmd.extend(
                ['vir_qual_filt.py', quality_summary_tsv, viral_posi_ctgs_fna, viral_posi_ctgs_fna, '\n']
            )
        else:
            cmd.extend(
                ['cp', viral_filt_ctgs_fna, ]
            )
        return cmd
    def Identify(self, cutoff = 1500, unrun = False):
        try:
            if int(self.allthreads) < 8:
                raise ValueError('The threads number must not be less than 8!!!')
        except ValueError as e:
            print(f'ERROR: {e}')
            exit(1)
        results = ''
        cutoff = str(cutoff)
        self.threads = str(self.threads)
        #vibrant
        cmd, wkdir = self.vibrant(cutoff)
        shell = f'{self.outdir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell, cmd)
        #deepvirfinder
        cmd, wkdir = self.deepvirfinder(cutoff)
        shell = f'{self.outdir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell, cmd)
        #genomad
        self.threads = str(int(self.threads))
        cmd, wkdir = self.genomad()
        shell = f'{self.outdir}/{self.name}_gm_ctg.sh'
        utils.printSH(shell, cmd)
        #VirSorter2
        self.threads = str(int(self.allthreads) - (int(self.threads) * 3))
        cmd, wkdir = self.virsorter(self.fasta, 0, cutoff)
        shell = f'{self.outdir}/{self.name}_vs2_ctg.sh'
        utils.printSH(shell, cmd)
        #multiple run
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(['multithreads.pl', self.outdir, 'ctg.sh 4\n'])
        shell = f'{self.outdir}/{self.name}_batch_identify_virus.sh'
        utils.printSH(shell, cmd)
        if not unrun: results = utils.execute(shell) 
        self.threads = str(int(self.allthreads))
        cmd = self.vFilter()
        shell = f'{self.outdir}/{self.name}_get_positive_virus.sh'
        utils.printSH(shell, cmd)
        if not unrun: results += utils.execute(shell)
        return results
