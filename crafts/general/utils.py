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
    SHPATH.write('#!/bin/bash\nset -e\n'+scripts)
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

def execute(sh_file:str):
    log_file=f'{sh_file}.log'
    error_file=f'{sh_file}.error'
    try:
        result=subprocess.run(['bash', sh_file],capture_output=True,text=True,check=True)
        with open(log_file,'w') as log_f:
            log_f.write(result.stdout)
        with open(error_file,'w') as error_f:
            error_f.write(result.stderr)
        return result.returncode
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a') as log_f:
            log_f.write(e.stdout)
        with open(error_file, 'a') as error_f:
            error_f.write(e.stderr)
        print(f"Script {sh_file} failed with error. Check {error_file} for details.")
        return e.returncode
    except Exception as e:
        print(f'Unexpected error: {e}')
        return 1

def run(cmd):
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

def is_file_exist(file_name):
    try:
        with open(file_name, 'r') as f:
            return 0
    except FileNotFoundError as e:
        print(f'FileNotFoundError: {e}')
        exit(1)
