import os
import re
import sys
import pandas as pd

def readConfig(config: str):
    '''
    Read the configure file and return a dictionary formated by "Key:value".
    '''
    confDict={}
    inconf=open(config)
    for line in inconf:
        if re.match(r'#|\s+',line): continue
        line=line.strip()
        [keys,values]=re.split(r'\s*=\s*',line,1)
        confDict[keys]=values
    inconf.close()
    #confDict['subProjectName']=re.split(r'_',confDict['ProjectName'])[1]
    return confDict

def readSampInfo(samp_info: str):
    ''' 
    '''
    sampDict={}
    insamp=open(samp_info)
    header=insamp.readline()
    header=header.strip().split('\t')
    id_idx=header.index('#SampleID')
    name_idx=header.index('SampleName')
    group_idx=header.index('GroupName')
    path_idx=header.index('DataPath')
    for line in insamp:
        items=line.strip().split('\t')
        sampDict[items[id_idx]]=[items[name_idx],items[group_idx],items[path_idx]]
    return sampDict

def prepInfo(config: str):
    confDict=readConfig(config)
    sampDict=readSampInfo(confDict['SampInfo'])
    return confDict,sampDict
    
