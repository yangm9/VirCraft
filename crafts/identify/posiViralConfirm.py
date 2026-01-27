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
        #self.threads = int(threads) // (self.BATCH_SIZE * 2) #self.BATCH_SIZE = 4, assiged in ..config.config
    def vFilter(self, min_len=2000, methods='gn', filt_mode='permissive'):
        score_tsv = f'{self.wkfile_dir}/all_viral_ctgs.score.tsv'
        score_filt_tsv = utils.insLable(score_tsv, filt_mode)
        viral_filt_ctg_list = f'{self.wkfile_dir}/viral_filt_ctg.list'
        viral_filt_ctgs_fna = f'{self.wkfile_dir}/viral_filt_ctg.fna'
        viral_posi_ctgs_fna = f'{self.wkfile_dir}/viral_positive_ctg.fna'
        venn_pdf = f'{self.stat_dir}/viral_contigs_multitools_comparison_venn.pdf'
        cmd = [utils.selectENV('VC-General')]
        tmp_cmd = ''
        #tmp_cmd, cat_dir = self.contig_annotation_tool(viral_filt_ctgs_fna)
        cmd.extend(
            ['merge_vctg_info.py', self.name, self.wkfile_dir, '-m', methods, '-f', filt_mode, # generate all_viral_ctgs.score.tsv file 
             '&& multitool_venn.py', '-i', score_tsv, '-t', methods, '-o', self.wkfile_dir, # plot Venn for all tools
             '&& cut -f 1', score_filt_tsv, "|sed '1d' >", viral_filt_ctg_list, 
             '&& extrSeqByName.pl', viral_filt_ctg_list, self.fasta, viral_filt_ctgs_fna, '\n\n']
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

    @staticmethod 
    def allocate_threads(enabled_tools, total_threads):
        BASE_WEIGHTS = {"vb": 1, "dvf": 1, "gn": 2, "vs2": 4}
        OVERRIDE_PAIRS = {frozenset({"gn", "vs2"}): {"gn": 2, "vs2": 6}}
        enabled_tools = set(enabled_tools)
        if len(enabled_tools) == 1:
            tool = next(iter(enabled_tools))
            return {tool: total_threads}
        allocation = {t: 1 for t in enabled_tools}
        remaining = total_threads - len(enabled_tools)
        if remaining < 0:
            raise ValueError("total_threads must be >= number of enabled tools")
        weights = {t: BASE_WEIGHTS[t] for t in enabled_tools}
        for pair, override in OVERRIDE_PAIRS.items():
            if pair.issubset(enabled_tools):
                for t in pair:
                    weights[t] = override[t]
        weight_sum = sum(weights.values())
        extra = {t: (remaining * w) // weight_sum for t, w in weights.items()}
        for t, v in extra.items():
            allocation[t] += v
        remainder = total_threads - sum(allocation.values())
        for t, _ in sorted(weights.items(), key=lambda x: x[1], reverse=True):
            if remainder == 0:
                break
            allocation[t] += 1
            remainder -= 1
        return allocation

    def Identify(self, min_len=2000, methods='gn', filt_mode='permissive', unrun=False):
        try:
            if int(self.threads) < 8:
                raise ValueError('The threads number must not be less than 8!!!')
        except ValueError as e:
            print(f'ERROR: {e}')
            exit(1)
        all_threads = self.threads
        enabled_tools = set(methods.split('-'))
        thread_map = self.allocate_threads(enabled_tools, int(self.threads))
        #vibrant
        if 'vb' in methods:
            self.threads = str(thread_map['vb'])
            cmd, wkdir = self.vibrant(min_len)
            shell = f'{self.shell_dir}/{self.name}_vb_ctg.sh'
            utils.printSH(shell, cmd)
        #deepvirfinder
        if 'dvf' in methods:
            self.threads = str(thread_map['dvf'])
            cmd, wkdir = self.deepvirfinder(min_len)
            shell = f'{self.shell_dir}/{self.name}_dvf_ctg.sh'
            utils.printSH(shell, cmd)
        #genomad
        if 'gn' in methods:
            self.threads = str(thread_map['gn'])
            self.threads = str(int(self.threads) * 2)
            cmd, wkdir = self.genomad()
            shell = f'{self.shell_dir}/{self.name}_gn_ctg.sh'
            utils.printSH(shell, cmd)
        #VirSorter2
        if 'vs2' in methods:
            self.threads = str(thread_map['vs2'])
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
        self.threads = str(all_threads) #8
        cmd = self.vFilter(min_len, methods, filt_mode)
        shell = f'{self.shell_dir}/{self.name}_get_positive_virus.sh'
        utils.printSH(shell, cmd)
        if not unrun: results += utils.execute(shell)
        return results
