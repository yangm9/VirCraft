import os
import sys
from ..general import utils
from ..votus.votus import VirRef

class VirHost(VirRef):
    '''
    Predict the relationship of virus and hosts.   
    '''
    def __init__(self,fasta='',hostsdir='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
        self.hostsdir=os.path.abspath(hostsdir)
    def magsTree(self):
        '''
        Classify the host MAGs by GTDBTK tools.
        '''
        wkdir=f'{self.outdir}/classify_wf'
        utils.mkdir(wkdir)
        cmd=['gtdbtk classify_wf','--genome_dir',self.hostsdir,
            '--out_dir',wkdir,'--pplacer_cpus',self.threads,
            '--extension fasta\n']
        return cmd,wkdir
    def virmatch(self,tredir):
        '''
        Predict the hosts for viral contigs using VirMatcher.
        '''
        wkdir=f'{self.outdir}/virmatcher'
        utils.mkdir(wkdir)
        cmd=['VirMatcher --preparer','--gtdbtk-out-dir',tredir,
            '--gtdbtk-in-dir',self.hostsdir,'-v',self.fasta,
            '-o',wkdir,'--threads',self.threads,'--python-aggregator']
        return cmd
    def PredHosts(self,unrun=False):
        cmd=[utils.selectENV('gtdbtk')]
        tmp_cmd,tredir=self.magsTree()
        cmd.extend(tmp_cmd)
        cmd.extend([self.envs])
        cmd.extend(self.virmatch(tredir))
        shell=f'{self.outdir}/{self.name}_hosts.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
