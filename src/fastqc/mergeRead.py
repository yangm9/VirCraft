import os
from ..config import conf
from ..process import cmdExec

def grpReads(config: str):
    global groups, confDict, sampDict
    groups, confDict, sampDict = conf.prepInfo(config)
    groups = []
    groupFq1Dict = {}
    groupFq2Dict = {}
    for samp in sampDict.keys():
        group = sampDict[samp][1]
        fastq1,fastq2 = sampDict[samp][2].split(',')
        groupFq1Dict.setdefault(group, [])
        groupFq2Dict.setdefault(group, [])
        groupFq1Dict[group].append(fastq1)
        groupFq2Dict[group].append(fastq2)
    return groupFq1Dict, groupFq2Dict

def catFastqCMD(fastqs: list, grp: str, n: int, outdir: str):
    merge_fq_cmd = ['cat']
    merge_fq_cmd.extend(fastqs)
    merge_dir = f'{outdir}/01.assembly/fastq'
    if not os.path.exists(merge_dir): os.makedirs(merge_dir)
    merge_fq_cmd.append(f'>{merge_dir}/{grp}_{n}.fq')
    return merge_fq_cmd

def MergeFastqByGroup(config: str,outdir: str):
    groupFq1Dict, groupFq2Dict = grpReads(config)
    for grp in groups:
        merge_fq1_cmd = catFastqCMD(groupFq1Dict[grp], grp, 1, outdir)
        merge_fq2_cmd = catFastqCMD(groupFq1Dict[grp], grp, 2, outdir)
        results1 = cmdExec.execute(merge_fq1_cmd)
        results2 = cmdExec.execute(merge_fq2_cmd)
    return results1, results2
