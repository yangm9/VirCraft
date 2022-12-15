import os
import re
import sys
from ..general import general,cmdExec

class VirCfg:
    '''
    Configure class.
    '''
    config=f'{sys.path[0]}/crafts/config/vircraft.cfg'
    def __init__(self):
        self.confDict=self.getCfg
    @property
    def getCfg(self):
        '''
        Read the configure file and return a dictionary formated by "Key:value".
        '''
        confDict={}
        inconf=open(self.config)
        for line in inconf:
            if re.match(r'#|\s+',line): continue
            line=line.strip()
            [keys,values]=re.split(r'\s*=\s*',line,1)
            confDict[keys]=values
        inconf.close()
        #confDict['subProjectName']=re.split(r'_',confDict['ProjectName'])[1]
        return confDict
    def readSampInfo(self,samp_info:str):
        '''
        Extract the groups list and the Sample information from "samp_info.xls" file. 
        Note: the format of sample information is "Sample: [GroupName,DataPath]".
        '''
        groups=[]
        sampDict={}
        insamp=open(samp_info)
        header=insamp.readline()
        header=header.strip().split('\t')
        name_idx=header.index('Sample')
        group_idx=header.index('Group')
        path_idx=header.index('DataPath')
        for line in insamp:
            items=line.strip().split('\t')
            sampDict[items[name_idx]]=[items[group_idx],items[path_idx]]
            groups.append(items[group_idx])
        groups=list(set(groups))
        insamp.close()
        return groups,sampDict

class Reads(VirCfg):
    '''
    FastQ processing class.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,fq1='',fq2='',outdir='',*args,**kwargs):
        super().__init__()
        self.fastqs=[fq1,fq2]
        self.basename_fq1=os.path.basename(self.fastqs[0])
        self.basename_fq2=os.path.basename(self.fastqs[1])
        self.outdir=os.path.abspath(outdir)
        self.samp=self.basename_fq1.replace('_1.fastq','')
        self.samp=self.samp.replace('_1.fq','')
        general.mkdir(self.outdir)

class Seq(VirCfg):
    '''
    Fasta processing class.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,fasta='',outdir='',*args,**kwargs):
        super().__init__()
        basename_fa=os.path.basename(fasta)
        self.name=os.path.splitext(basename_fa)[0]
        self.fasta=os.path.abspath(fasta)
        self.outdir=os.path.abspath(outdir)
        general.mkdir(self.outdir)
    @property
    def mkBwaIdx(self):
        "Make bwa index for votus."
        cmd=[self.envs]
        idx=f'{self.outdir}/{self.name}BWAIDX'
        cmd.extend(['bwa index -a bwtsw',self.fasta,'-p',idx,'\n'])
        shell=f'{self.outdir}/{self.name}_bwaidx.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return idx,results
