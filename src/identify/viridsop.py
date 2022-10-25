
from ..config import setVari,conf
from ..process import cmdExec

def findVir(grp: str, outdir: str):
    envs = setVari.selectENV('VirCraft')
    virsorter_cmd = [envs]
    identify_dir = f'{outdir}/02.identify/{grp}'
    scaffolds = f'{outdir}/01.assembly/{grp}/scaffolds.fasta'
    if not os.path.exists(identify_dir): os.makedirs(identify_dir)
    virsorter_cmd.extend(['virsorter', 'run', '--keep-original-seq', '-i', scaffolds, '-d ', '--include-groups dsDNAphage,ssDNA --min-length 0 --min-score 0.5 -j 28 all', '-o', spades_dir])

