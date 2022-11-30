import os
from ..general import cmdExec,general
from ..config.cfgInfo import VirCfg

class Reads(VirCfg):
    '''
    FastQ processing class.
    '''
    def __init__(self,config,outdir):
        VirCfg.__init__(self,config)
        self.outdir=os.path.abspath(outdir)
        self.wkdir=f'{self.outdir}/00.data'
        general.mkdir(self.wkdir)
        self.fq_dir=f'{self.wkdir}/merged_fq'
        general.mkdir(self.fq_dir)
        self.grpFq1Dict,self.grpFq2Dict=self.grpReads
    @property
    def grpReads(self):
        groupFq1Dict={}
        groupFq2Dict={}
        for samp in self.sampDict.keys():
            group=self.sampDict[samp][1]
            fastq1,fastq2=self.sampDict[samp][2].split(',')
            groupFq1Dict.setdefault(group,[])
            groupFq2Dict.setdefault(group,[])
            groupFq1Dict[group].append(fastq1)
            groupFq2Dict[group].append(fastq2)
        return groupFq1Dict,groupFq2Dict
    def catFqCMD(self,fastqs:list,grp:str,n:int):
        merge_fq_cmd=['cat']
        merge_fq_cmd.extend(fastqs)
        merge_fq_cmd.append(f'>{self.fq_dir}/{grp}_{n}.fq\n')
        return merge_fq_cmd
    @property
    def MergeFastqByGroup(self):
        results1,results2='',''
        for grp in self.groups:
            merge_fq1_cmd=self.catFqCMD(self.grpFq1Dict[grp],grp,1)
            merge_fq2_cmd=self.catFqCMD(self.grpFq2Dict[grp],grp,2)
            merge_fq_sh=f'{self.wkdir}/{grp}_merge_fq.sh'
            general.printSH(merge_fq_sh,merge_fq1_cmd+merge_fq2_cmd)
            results1 += cmdExec.execute(merge_fq1_cmd)
            results2 += cmdExec.execute(merge_fq2_cmd)
        return results1,results2
