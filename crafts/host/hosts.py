import os
import sys
from ..general import utils
from ..votus.uniqVirCtg import VirRef

class VirHost(VirRef):
    '''
    Predict the relationship of virus and hosts.   
    '''
    def __init__(self, fasta=None, hostsdir=None, outdir=None, threads=8):
        super().__init__(fasta, outdir, threads)
        self.hostsdir = os.path.abspath(hostsdir)
    def magsTree(self):
        '''
        Classify the host MAGs by GTDBTK tools.
        '''
        wkdir = f'{self.wkfile_dir}/classify_wf'
        utils.mkdir(wkdir)
        cmd = [utils.selectENV('VC-GTDBTk')]
        cmd.extend(
            ['gtdbtk classify_wf', '--genome_dir', self.hostsdir, '--out_dir', wkdir, '--pplacer_cpus', self.threads, '--extension fasta --skip_ani_screen\n']
        )
        return cmd, wkdir
    def virmatch(self, tredir):
        '''
        Predict the hosts for viral contigs using VirMatcher.
        '''
        wkdir = f'{self.wkfile_dir}/virmatcher'
        utils.mkdir(wkdir)
        cmd = [utils.selectENV('VC-VHMatcher')]
        cmd.extend(
            ['VirMatcher --preparer', '--gtdbtk-out-dir', tredir, '--gtdbtk-in-dir', self.hostsdir, '-v', self.fasta, '-o', wkdir, '--threads', self.threads, '--python-aggregator\n']
        )
        return cmd, wkdir
    def iphop(self, tredir: str):
        '''
        Predict the hosts for viral contigs using VirMatcher.
        '''
        wkdir = f'{self.wkfile_dir}/iphop'
        utils.mkdir(wkdir)
        cmd = [utils.selectENV('VC-iPHoP')]
        cmd.extend(
            ['iphop add_to_db', '--fna_dir', self.hostsdir, '--gtdb_dir', tredir, '--db_dir', self.confDict['iPHoPDB'], '--num_threads', self.threads, '\n'],
        )
        return cmd, wkdir
    def virTaxa(self, taxa_anno=None):
        '''
        Add the viral taxonomic annotation for host pridiction results.
        '''
        wkdir = f'{self.wkfile_dir}/virmatcher'
        m_taxa_anno = f'{self.wkfile_dir}/all_votu.taxa.txt'
        vh_pred = f'{wkdir}/VirMatcher_Summary_Predictions.tsv'
        taxa_anno = os.path.abspath(taxa_anno)
        vh_vtaxa = f'{self.wkfile_dir}/VirMatcher_Summary_Predictions.vtaxa.tsv'
        if not taxa_anno:
            return ['ln -s', vh_pred, vh_vtaxa, '\n']
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ["sed -i '1s/ /_/g'", vh_pred, "&& sed '1s/Sequence_ID/Original_Viral_population/'", taxa_anno, '>', m_taxa_anno, '\n',
             'linkTab.py', vh_pred, m_taxa_anno, 'left Original_Viral_population', vh_vtaxa,'\n']
        )
        return cmd
    def hostTaxa(self, tredir):
        tredir = os.path.abspath(tredir)
        gtdbtk_arc = f'{tredir}/gtdbtk.ar53.summary.tsv'
        gtdbtk_bac = f'{tredir}/gtdbtk.bac120.summary.tsv'
        gtdbtk_bac_tmp = gtdbtk_bac + '.tmp'
        gtdbtk_hosts = f'{self.wkfile_dir}/gtdbtk.hosts.taxa.tsv'
        vh_vtaxa = f'{self.wkfile_dir}/VirMatcher_Summary_Predictions.vtaxa.tsv'
        vh_vhtaxa = vh_vtaxa.replace('vtaxa', 'vhtaxa')
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ["sed '1d'", gtdbtk_bac, '>', gtdbtk_bac_tmp, '\n',
             'cat', gtdbtk_arc, gtdbtk_bac_tmp, "|cut -f 1,2|sed '1s/user_genome/Original_Host/' >", gtdbtk_hosts, '\n',
             'linkTab.py', vh_vtaxa, gtdbtk_hosts, 'left Original_Host', vh_vhtaxa, '\n'])
        return cmd
    def PredHosts(self, gtdbtk=None, taxa_anno=None, unrun=False):
        cmd = []
        if not gtdbtk:
            tmp_cmd,gtdbtk=self.magsTree()
            cmd.extend(tmp_cmd)
        cmd.extend(self.virmatch(gtdbtk))
        cmd.extend(self.virTaxa(taxa_anno))
        cmd.extend(self.hostTaxa(gtdbtk))
        shell = f'{self.shell_dir}/{self.name}_hosts.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
