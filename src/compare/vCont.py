from ..config import setVari, conf
from ..process import cmdExec, general

def compSeq(config: str, outdir: str):
    groups, confDict, sampDict = conf.prepInfo(config)
    envs = setVari.selectENV('vContact2')
    comp_cmd = [envs]
    comp_dir = f'{outdir}/07.comparison/'
    prodigal_dir = '{comp_dir}/0.prodigal'
    vcont_dir = '{comp_dir}/1.vContact2'
    general.mkdir(comp_dir)
    general.mkdir(prodigal_dir)
    general.mkdir(vcont_dir)
    orfs_faa = f'{outdir}/06.func_annot/0.prodigal/merged_virus_positive_nodup.faa'
    merged_faa = f'{prodigal_dir}/all.merged.fa'
    prot_info_csv = f'{prodigal_dir}/all.merged.csv'
    vcontact2DB = "'ProkaryoticViralRefSeq94-Merged'"
    comp_cmd.expend(
        ['cat', orfs_faa, confDict['OtherDataSetFAA'], '>', merged_faa, '\n',
        'virus2csv.py', merged_faa, '>', prot_info, '\n',
        'vcontact2 --rel-mode Diamond --pcs-mode MCL --vcs-mode ClusterONE -t 28', '--raw-proteins', merged_faa, '--proteins-fp', prot_info_csv, '--db', vcontactDB, '--output-dir', vcont_dir, '\n']
    )
    comp_sh = f'{comp_dir}/compare.sh'
    general.printSH(comp_sh, comp_cmd)
    results = cmdExec.execute(comp_cmd)
    return results
