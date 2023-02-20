import os
import sys
from ..general import cmdExec,general
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
        general.mkdir(wkdir)
        cmd=['gtdbtk classify_wf','--genome_dir',self.hostsdir,
            '--out_dir',wkdir,'--pplacer_cpus',self.threads,
            '--extension fasta\n']
        return cmd,wkdir
    def virmatch(self,tredir):
        '''
        Predict the hosts for viral contigs using VirMatcher.
        '''
        wkdir=f'{self.outdir}/virmatcher'
        general.mkdir(wkdir)
        cmd=['VirMatcher --preparer','--gtdbtk-out-dir',tredir,
            '--gtdbtk-in-dir',self.hostsdir,'-v',self.fasta,
            '-o',wkdir,'--threads',self.threads,'--python-aggregator']
        return cmd
    def PredHosts(self):
        cmd=[general.selectENV('gtdbtk')]
        tmp_cmd,tredir=self.magsTree()
        cmd.extend(tmp_cmd)
        cmd.extend([self.envs])
        cmd.extend(self.virmatch(tredir))
        shell=f'{self.outdir}/{self.name}_hosts.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
