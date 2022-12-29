import os
from ..general import cmdExec,general
from ..config.config import ORF

class Dataset:
    def __init__(self,fa_list='',outdir='',threads=8):
        self.fa_list=os.path.abspath(fa_list)

class EnviComp(ORF):
    '''
    Compare the protein sequence from current to other environments.
    '''
    def __init__(self,orfs='',outdir='',threads=8):
        super().__init__(orfs,outdir,threads)
        self.envs=general.selectENV('vContact2')
    def compSeq(self):
        cmd=[self.envs]
        wkdir=f'{self.outdir}/vContact2'
        general.mkdir(wkdir)
        #OtherORFs=self.confDict['OtherDataSetFAA']
        vContDB="'ProkaryoticViralRefSeq94-Merged'"
        orfs_info=f'{wkdir}/all_merged_orfs.csv'
        comp_cmd.expend(
            ['virus2csv.py',self.orfs,'>',orfs_info,'\n',
            'vcontact2 --rel-mode Diamond --pcs-mode MCL --vcs-mode ClusterONE',
            '-t',self.threads,'--raw-proteins',merged_orfs,
            '--proteins-fp',orfs_info,'--db',vContDB,
            '--output-dir',self.outdir,'\n']
        )
        shell=f'{comp_dir}/v_compare.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
