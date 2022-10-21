import os
import re
import sys

def readConfig(config):
    'Read the configure file and return a dictionary formated by "Key:value".'
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
