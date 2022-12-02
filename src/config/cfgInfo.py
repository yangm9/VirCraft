import os
import re
import sys

class VirCfg:
    '''
    Configure class.
    '''
    config=f'{sys.path[0]}/src/config/vircraft.cfg'
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
