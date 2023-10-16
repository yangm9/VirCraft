import os
import sys
import importlib
import subprocess
import logging
from datetime import datetime 

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
    conda_path_dirs=conda_path.split('/')
    mc3_dir_idx=conda_path_dirs.index('miniconda3')
    condash_path='/'.join(conda_path_dirs[0:mc3_dir_idx+1])
    condash_path+='/etc/profile.d/conda.sh'
    envs=''
    if os.path.exists(condash_path):
        envs=f'source "{condash_path}"\nconda activate && conda activate {env}\n'
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

def run(cmd:list,sh_file:str,silent=False):
    cmd_txt = ' '.join(cmd).replace('\n ','\n')
    cmd_txt = cmd_txt.replace(' \n ','\n')
    print(f'Running command:\n{cmd_txt}') 
    if not silent:
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, shell=False, check=True)
    if result.returncode != 0:
        print(f'Error running command: {cmd_txt}. The error message was: \n{result.stderr}')
        exit(1)
    return result

def execute(cmd):
    cmd_txt=' '.join(cmd).replace('\n ','\n')
    cmd_txt=cmd_txt.replace(' \n ','\n')
    print(f'Running command:\n{cmd_txt}')
    results=str(os.system(cmd_txt))
    return results

def show_cmd(cmd):
    cmd_txt=' '.join(cmd).replace('\n ','\n')
    cmd_txt=cmd_txt.replace(' \n ','\n')
    print(f'Pending command:\n{cmd_txt}')
    return 0

def install_module(module_name):
    try:
        module=importlib.import_module(module_name)
        print(f"Module '{module_name}' is already installed.")
        return module
    except ImportError:
        print(f"Module '{module_name}' is not installed. Attempting to install...")
        try:
            subprocess.check_call(['pip3', 'install', module_name])
            print(f"Module '{module_name}' installed successfully.")
        except subprocess.CalledProcessError:
            print(f"Failed to install module '{module_name}'.")
            return None
        try:
            module = importlib.import_module(module_name)
            print(f"Module '{module_name}' is now installed and imported.")
            return module
        except ImportError:
            print(f"Module '{module_name}' was installed but cannot be imported.")
            return None

try:
    import warnings
    from conda.base.context import context
    warnings.filterwarnings("ignore")
except ImportError as e:
    print(f'''ImportError: {e}
Please activate conda base environment using `conda activate` command...''')
    exit(1)
#    try:
#        subprocess.run(['conda','activate','base'],shell=True)
#        from conda.base.context import context
#    except Exception as activation_error:
#        print("Failed to activate Conda environment.")
#        print(str(activation_error))
#    #install_module('conda')

def get_conda_env_dir(env_name):
    envs_dir=os.path.join(context.root_prefix,'envs')
    env_dir=os.path.join(envs_dir,env_name)
    return env_dir
