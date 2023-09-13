import sys
from datetime import datetime
from ..general import utils
from . import URL

class ENV:
    '''
    Install the environments of VirCraft
    '''
    CONDAENVDIR=f'{sys.path[0]}/crafts/setup/conda_env_yaml'
    ENVLIST=['reads_qc','assembly','viral-id-sop','vibrant','deepvirfinder',
            'vcontact','kofamscan','gtdbtk','kofamscan','vircraft']
    def __init__(self,outdir='',threads=8):
        self.outdir=outdir
        self.threads=str(threads)
    def setup_env(self,name,in_wall=False):
        net=''
        if in_wall: net='_cn'
        name=name+net+'.yaml'
        cmd=['conda env create','-f',name,'\n']
        return cmd
    def Install(self,in_wall=False):
        for env in ENVLIST:
            cmd=self.setup_env(env,in_wall)
            shell=f'{self.outdir}/{env}_install.sh'
            utils.printSH(shell,cmd)
            results=''
            if not unrun: results+=utils.execute(cmd)
        return results
