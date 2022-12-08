import os
from ..general import cmdExec,general
from ..func_annot.geneAnnot import GeneFunc

class Dataset:
    def __init__(self,fa_list='',outdir='',threads=8):
        self.fa_list=os.path.abspath(fa_list)

        

class EnviComp(GeneFunc):
    '''
    Compare the protein sequence from current to other environments.
    '''
    def __init__(self,orfs='',outdir='',threads=8):
        super().__init__(orfs,outdir,threads)
    def compSeq(self,prot_dataset):
        cmd=[self.envs]
        general.mkdir(prodi_dir)
        vcont_dir='{self.wkdir}/1.vContact2'
        general.mkdir(vcont_dir)
        vContDB="'ProkaryoticViralRefSeq94-Merged'"
        comp_cmd.expend(
            ['cat',orfs_faa,confDict['OtherDataSetFAA'],'>',merged_faa,'\n',
            'virus2csv.py',merged_faa,'>',prot_info,'\n',
            'vcontact2 --rel-mode Diamond --pcs-mode MCL --vcs-mode ClusterONE,',
            '-t',self.threads,'--raw-proteins',merged_faa,
            '--proteins-fp',prot_info,'--db',vContDB,
            '--output-dir',vcont_dir,'\n']
        )
        shell=f'{comp_dir}/v_compare.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
