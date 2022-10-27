#!/usr/bin/env python3

from os import path
from ..config import setVari, conf
from ..process import cmdExec, general

def seqCluster(fasta: str, group: str, outdir: str):
    '''
    Cluster the sequence and remove redundancy for FastA file.
    '''
    envs = setVari.selectENV('VirCraft')
    cdhit_cmd = [envs]
    fa_prefix=path.splitext(path.basename(fasta))[0]
    nodup_fa = f'{outdir}/{group}_{fa_prefix}_nodup.fa'
    cdhit_cmd.extend(
        ['cd-hit-est', '-i', fasta, '-o', nodup_fa, 
        '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000 -T 28\n']
    )
    cdhit_sh = f'{outdir}/{group}_cdhit.sh'
    general.printSH(cdhit_sh, cdhit_cmd)
    results = cmdExec.execute(cdhit_cmd)
    return results

def RmDup(config: str, outdir: str):
    groups, confDict, sampDict = conf.prepInfo(config)
    results = ''
    wkdir = f'{outdir}/03.vOTUs'
    general.mkdir(wkdir)
    merge_fa_cmd = ['cat']
    for grp in groups:
        fasta = f'{outdir}/02.identify/{grp}/curation/viruses_positive.fna'
        merge_fa_cmd.append(fasta)
        grpdir = f'{wkdir}/{grp}'
        general.mkdir(grpdir)
        results += f'{grp}: \n'
        results += seqCluster(fasta, grp, wkdir)
    merged_fasta = f'{wkdir}/merged_virus_positive_nodup.fa'
    merge_fa_cmd.append(f'>{merged_fasta}\n')
    merge_fa_sh = f'{wkdir}/merge_fa.sh'
    general.printSH(merge_fa_sh,merge_fa_cmd)
    results += cmdExec.execute(merge_fa_cmd)
    results += seqCluster(merged_fasta, 'merged', wkdir)
    return results
