import sys
import subprocess
from datetime import datetime
from ..general import utils
from . import URL

class ENV:
    '''
    Install the environments of VirCraft
    '''
    CONDAENVDIR=f'{sys.path[0]}/crafts/setup/conda_env_yaml'
    ENVDICT={'reads_qc':'VC-ReadsQC','assembly':'VC-Assembly',
             'virsorter2':'VC-VirSorter2','vibrant':'VC-VIBRANT',
             'deepvirfinder':'VC-DeepVirFinder','checkv':'VC-CheckV',
             'dramv':'VC-DRAMv','vcontact':'VC-vContact2',
             'gtdbtk':'VC-GTDBTk','quantify':'VC-Quantify',
             'vhmatcher':'VC-VHMatcher','general':'VC-General'}
    def __init__(self,outdir='',threads=8):
        self.outdir=outdir
        self.threads=str(threads)
        utils.mkdir(self.outdir)
    def is_conda_env(self,env_name):
        try:
            # 使用conda命令列出已安装的环境
            result=subprocess.run(['conda','env','list'],capture_output=True,text=True,check=True)
            output_lines=result.stdout.splitlines()
            # 检查每行输出是否包含指定的环境名称
            env_exist=any(env_name in line for line in output_lines)
            if env_exist:
                return 1
            else:
                return 0
        except subprocess.CalledProcessError as e:
            # 如果无法执行conda命令，则返回-1表示出错
            print(f'ERROR: {e}')
            exit(1)
            return -1
    def setup_env(self,name,in_wall=False):
        net=''
        if in_wall: net='_cn'
        env_yaml=f'{self.CONDAENVDIR}/{name}{net}.yaml'
        cmd=['mamba env create','-f',env_yaml,'\n']
        return cmd
    def Install(self,in_wall=False,unrun=False):
        for env in self.ENVDICT.keys():
            is_env=self.is_conda_env(self.ENVDICT[env])
            if is_env:
                print(f'{self.ENVDICT[env]} is already installed, skipping!')
                continue
            cmd=self.setup_env(env,in_wall)
            if env=='':
                wishdir=''
                cmd.extend(
                    ['git clone',URL.WISH_URL,'cd',wishdir,
                    '&& cmake . && make && chmod +x WIsH && cp WIsH',]
                )
            shell=f'{self.outdir}/{env}_install.sh'
            utils.printSH(shell,cmd)
            results=''
            if not unrun: results+=utils.execute(cmd)
        return results
