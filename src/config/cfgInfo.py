import os
import re
import sys
import pandas as pd

class VirCfg:
    '''
    Configure class.
    '''
    def __init__(self, config):
        self.config = os.path.abspath(config)
        self.confDict = self.readConfig
        self.samp_info = os.path.abspath(self.confDict['SampInfo'])
        self.groups, self.sampDict = self.readSampInfo
    @property
    def readConfig(self):
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
    @property
    def readSampInfo(self):
        '''
        Extract the groups list and the Sample information from "samp_info.xls" file. 
        Note: the format of sample information is "Sample ID : [SampleName,GroupName,DataPath]".
        '''
        groups=[]
        sampDict={}
        insamp=open(self.samp_info)
        header=insamp.readline()
        header=header.strip().split('\t')
        id_idx=header.index('SampleID')
        name_idx=header.index('SampleName')
        group_idx=header.index('GroupName')
        path_idx=header.index('DataPath')
        for line in insamp:
            items=line.strip().split('\t')
            sampDict[items[id_idx]]=[items[name_idx],items[group_idx],items[path_idx]]
            groups.append(items[group_idx])
        groups=list(set(groups))
        insamp.close()
        return groups,sampDict