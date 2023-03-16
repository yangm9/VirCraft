import os
import sys

def mkdir(name:str):
    if not os.path.exists(name):
        os.makedirs(name)
    return 0

def printSH(sh_path:str,command:list):
    SHPATH=open(sh_path,'w')
    scripts=' '.join(command).replace('\n ','\n')
    scripts=scripts.replace(' \n ','\n')
    SHPATH.write(scripts)
    SHPATH.close()
    return 0

def insLable(file_name:str,label:str):
    "Insert a label before the extension for a file name."
    prefix,extension=os.path.splitext(file_name)
    return f'{prefix}.{label}{extension}'

def selectENV(env:str):
    bin_dir=os.path.abspath(sys.path[0]+'/bin')
    conda_path=isInstalled('conda')
    condash_path='/'.join(conda_path.split('/')[0:-2])
    condash_path+='/etc/profile.d/conda.sh'
    envs=''
    if os.path.exists(condash_path):
        envs=f'. "{condash_path}"\nconda activate {env}\n'
    else:
        conda_bin=os.path.dirname('conda_path')
        envs=f'export PATH="{conda_bin}:$PATH"\n'
    envs+=f'export PATH="{bin_dir}:$PATH"\n'
    return envs
def isInstalled(name:str):
    """
    Check whether `name` is on PATH and marked as executable.
    """
    from shutil import which
    return which(name)
