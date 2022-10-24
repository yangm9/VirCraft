from ..process import cmdExec

def grpReads(config: str):
    confDict,sampDict=conf.prepInfo(config)
    groups=[]
    groupFq1Dict={}
    groupFq2Dict={}
    for samp in sampDict.keys():
        group=sampDict[samp][1]
        groups.append(group)
        fastq1,fastq2=sampDict[samp][2].split(',')
        groupFq1Dict.setdefault(group,[])
        groupFq2Dict.setdefault(group,[])
        groupFq1Dict[group].append(fastq1)
        groupFq2Dict[group].append(fastq2)
    return groups,groupFq1Dict,groupFq2Dict

def catFastqCMD(fastqs: list,n: int,outdir: str):
    merge_fq_cmd=['cat']
    merge_fq_cmd.extend(fastqs)
    merge_fq_cmd.append(f'>{outdir}/01.assembly/0.fastq/{grp}_{n}.fq}')
    return merge_fq_cmd

def MergeFastqByGroup(config: str,outdir: str):
    groups,groupFq1Dict,groupFq2Dict=grpReads(config)
    for grp in groups:
        merge_fq1_cmd=catFastqCMD(groupFq1Dict[grp],1,outdir)
        merge_fq1_cmd=catFastqCMD(groupFq1Dict[grp],2,outdir)
        results1=cmdExec.execute(merge_fq1_cmd)
        results2=cmdExec.execute(merge_fq2_cmd)
    return results1,results2
