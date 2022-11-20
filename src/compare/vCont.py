from ..general import cmdExec,general
from ..func_annot.geneAnnot import GeneFunc
from ..votus.deRep import VirRef 

class EnviComp(VirRef):
    '''
    '''
    def __init__(self,config,outdir):
        VirRef.__init__(self,config,outdir)
        self.datadir=self.wkdir
        self.wkdir=f'{self.outdir}/07.compare'
        general.mkdir(self.wkdir)
    def compSeq(self):
        cmd=[self.envs]
        prodi_dir='{self.wkdir}/0.prodigal'
        general.mkdir(prodi_dir)
        vcont_dir='{self.wkdir}/1.vContact2'
        general.mkdir(vcont_dir)
        merged_faa=f'{prodi_dir}/all_merged.faa'
        prot_info=f'{prodi_dir}/all_merged.csv'
        vcontact2DB="'ProkaryoticViralRefSeq94-Merged'"
        comp_cmd.expend(
            ['cat',orfs_faa,confDict['OtherDataSetFAA'],'>',merged_faa,'\n',
            'virus2csv.py',merged_faa,'>',prot_info,'\n',
            'vcontact2 --rel-mode Diamond --pcs-mode MCL --vcs-mode ClusterONE -t 28','--raw-proteins',merged_faa,'--proteins-fp',prot_info,'--db',vcontactDB,'--output-dir',vcont_dir,'\n']
        )
        shell=f'{comp_dir}/v_compare.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
