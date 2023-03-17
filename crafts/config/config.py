import os
import re
import sys
from ..general import utils

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
        Extract the groups list and the Sample information from "samp_info.xls" file, with the format: "Sample\tGroup\tDataPath\n". 
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
