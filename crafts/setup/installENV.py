import os
import sys
import subprocess
from datetime import datetime
from ..general import utils
from . import URL

class ENV:
    '''
    Install the environments of VirCraft
    '''

    CONDAENVDIR = f'{sys.path[0]}/crafts/setup/conda_env_yaml'
    ENVDICT = {'reads_qc': 'VC-ReadsQC', 'assemble': 'VC-Assemble', 'virsorter2': 'VC-VirSorter2', 'vibrant': 'VC-VIBRANT', 
               'deepvirfinder': 'VC-DeepVirFinder', 'genomad': 'VC-geNomad', 'checkv': 'VC-CheckV', 'vrhyme': 'VC-vRhyme',
               'dramv': 'VC-DRAMv', 'vcontact': 'VC-vContact2', 'gtdbtk': 'VC-GTDBTk', 'iphop': 'VC-iPHoP', 'quantify': 'VC-Quantify', 
               'vhmatcher': 'VC-VHMatcher', 'general': 'VC-General'}
    def __init__(self, outdir=None, threads=8):
        self.outdir = os.path.abspath(outdir)
        self.threads = str(threads)
        utils.mkdir(self.outdir)
    def is_conda_env(self, env_name):
        try:
            result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True, check=True)
            output_lines = result.stdout.splitlines()
            env_exist = any(env_name in line for line in output_lines)
            if env_exist:
                return 1
            else:
                return 0
        except subprocess.CalledProcessError as e:
            print(f'ERROR: {e}')
            exit(1)
            return -1
    def setup_env(self, name, in_wall=False):
        net = ''
        if in_wall: net = '_cn'
        env_yaml = f'{self.CONDAENVDIR}/{name}{net}.yaml'
        cmd = [utils.selectENV('base')]
        cmd.extend(['mamba env create', '-f', env_yaml, '-y\n'])
        if name == 'vhmatcher':
            tmp_wishdir = f'{self.outdir}/WIsH'
            vhmatcher_bin_dir = utils.get_conda_env_dir(self.ENVDICT[name])
            vhmatcher_bin_dir += '/bin'
            cmd.extend(
                ['mkdir', tmp_wishdir, '&& git clone', URL.WISH_URL, tmp_wishdir,
                 '&& cd', tmp_wishdir, '&& cmake . && make && chmod +x WIsH',
                 '&& cp WIsH', vhmatcher_bin_dir, '\n']
            )
            tmp_virmatcherdir = f'{self.outdir}/VirMatcher'
            cmd.extend(
                ['mkdir', tmp_virmatcherdir, '&& git clone', URL.VIRMATCHER_URL, tmp_virmatcherdir,
                 '&& cd', tmp_virmatcherdir, "&& sed -i 's/4.2_5/4/' setup.py", "&& sed -i 's/ar122/ar53/' bin/VirMatcher",
                 "&& sed -i -E 's/^indexes = set\(\).union\(\*indices_to_use\)/indexes = list(set().union(*indices_to_use))/' bin/ResultsAggregator.py",
                 "&& sed -i -E 's/^columns = set\(\).union\(\*columns_to_use\)/columns = list(set().union(*columns_to_use))/' bin/ResultsAggregator.py",
                 '&& conda run -n', self.ENVDICT[name], 'pip install . --no-deps\n']
            )
        elif name == 'general':
            tmp_catdir = f'{self.outdir}/CAT'
            cat_pack_files = f'{tmp_catdir}/CAT_pack/*'
            general_bin_dir = utils.get_conda_env_dir(self.ENVDICT[name])
            general_bin_dir += '/bin'
            cmd.extend(
                ['git clone', URL.CAT_URL, tmp_catdir,
                 '&& cp', cat_pack_files, general_bin_dir, '\n']
            )
        elif name == 'vrhyme':
            tmp_vrhymedir = f'{self.outdir}/vRhyme'
            model_gz = f'{tmp_vrhymedir}/models/vRhyme_machine_model_ET.sav.gz'
            cmd.extend(
                ['git clone', URL.VRHYME_URL, tmp_vrhymedir,
                 '&& cd', tmp_vrhymedir, '&& gunzip', model_gz,
                 '&& conda run -n', self.ENVDICT[name], 'pip install .\n']
            )
        return cmd
    def Install(self, in_wall=False, unrun=False):
        results = ''
        for env in self.ENVDICT.keys():
            is_env = self.is_conda_env(self.ENVDICT[env])
            cmd = self.setup_env(env, in_wall)
            shell = f'{self.outdir}/{env}_install.sh'
            utils.printSH(shell, cmd)
            if is_env:
                print(f'{self.ENVDICT[env]} installed, skipped!')
                continue
            if not unrun:
                print(f'{self.ENVDICT[env]} pending...')
                results += str(utils.execute(shell))
        return results
