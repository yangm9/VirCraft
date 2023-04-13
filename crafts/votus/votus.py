from os import path
from ..general import utils
from ..identify.viridsop import VirScan

class VirRef(VirScan):
    '''
    Generate vOTUs
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
    def cluster(self):
        '''
        Cluster the sequence and remove redundancy for FastA file.
        '''
        votus=f'{self.outdir}/{self.name}_votus.fa'
        cmd=['cd-hit-est','-i',self.fasta,'-o',votus,'-T',self.threads,
            '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000\n']
        return cmd,votus
    def RmDup(self,cutoff=5000,unrun=False):
        cmd=[self.envs]
        tmp_cmd,votus=self.cluster()
        cmd.extend(tmp_cmd)
        tmp_cmd,tmp_fa=self.checkv(votus)
        cmd.extend(tmp_cmd)
        tmp_cmd=self.statFA(cutoff)
        cmd.extend(tmp_cmd)
        shell=f'{self.outdir}/{self.name}_votu.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
