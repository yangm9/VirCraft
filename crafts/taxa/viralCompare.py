import os
from ..general import utils
from ..data.bioseq import ORF

class Dataset:
    def __init__(self, fa_list=None, outdir=None, threads=8):
        self.fa_list = os.path.abspath(fa_list)

class EnviComp(ORF):
    '''
    Compare the protein sequence from current to other environments.
    '''
    def __init__(self, fasta=None, outdir=None, threads=8, orfprefix=None):
        super().__init__(fasta, outdir, threads)
        self.threads = threads
        self.orfprefix = orfprefix
        #the first element of array generated by splitting orfs ID.
    def prodigal(self):
        return 0
    def vcontact2(self, vcont_db='None'):
        cmd = [utils.selectENV('VC-vContact2')]
        wkdir = f'{self.wkfile_dir}/vcontact2'
        utils.mkdir(wkdir)
        orfs_info = f'{self.wkfile_dir}/orf_to_contig.csv'
        cmd.extend(
            ['connect_orf_to_contig.py', self.fasta, '>', orfs_info, '\n',
             'vcontact2 --rel-mode Diamond --pcs-mode MCL --vcs-mode ClusterONE', '-t', self.threads, '--raw-proteins', self.fasta, '--proteins-fp', orfs_info, '--db', vcont_db, '--output-dir', wkdir, '\n']
        )
        if self.orfprefix:
            ntw = f'{wkdir}/c1.ntw'
            filt_ntw = f'{wkdir}/c1.filt.ntw'
            cmd.extend(
                ['grep', '"', self.orfprefix, '"', ntw, '>', filt_ntw, '\n']
            )
        cmd.extend(
            []
        )
        return cmd
    def CompSeq(self, unrun=False):
        cmd = self.vcontact2('')
        shell = f'{self.shell_dir}/viral_datasets_compare.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
