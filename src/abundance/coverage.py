from ..config import setVari,conf
from ..process import cmdExec, general

def calcTPM(samp: str, sort_bam: str):
    envs = setVari.selectENV('VirCraft')
    tpm_cmd = [envs]
    tpm = f'{wkdir}/{samp}.tpm'
    tpm_cmd.extend(
        ['coverm contig', '-b', sort_bam, 
        '-t 20 --min-read-aligned-length 50 --min-read-percent-identity 0.95 --proper-pairs-only -m tpm', 
        '>', tpm]
    )

    tpm_sh = f'{wkdir}/{samp}_tpm.sh'
    general.printSH(bwa_idx_sh, tpm_cmd)
    results = cmdExec.execute(tpm_cmd)
    return results

def TPMBySamp(config: str, outdir: str):
    global wkdir
    groups, confDict, sampDict = conf.prepInfo(config)
    wkdir = f'{outdir}/05.abundance/2.coverm'
    bam_dir = f'{outdir}/05.abundance/1.bwa'
    general.mkdir(wkdir)
    for samp in sampDict.keys():
        sort_bam = f'{bam_dir}/{samp}.sort.bam'
        results += calcTPM(samp, sort_bam)
    return results
