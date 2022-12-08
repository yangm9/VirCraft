import sys
from ..process import cmdExec,general
from ..votus.deRep import VirRef

class VirHost(VirRef):
    '''
    Predict the relationship of virus and hosts.   
    '''
    def __init__(self,config,outdir):
        VirRef.__init__(self,config,outdir)
        self.datadir=self.wkdir
        self.wkdir=f'{self.outdir}/08.hosts'
        self.tredir=f'{self.wkdir}/classify_wf'
        self.host_mags=self.confDict['HostMAGs']
        general.mkdir(self.wkdir)
    def magsTree(self):
        '''
        Classify the host MAGs by GTDBTK tools.
        '''
        cmd=[self.envs]
        cmd.extend(
            ['gtdbtk classify_wf','--genome_dir',self.host_mags,
            '--out_dir',self.tredir,'--extension fasta --pplacer_cpus 32\n']
        )
        shell=f'{self.wkdir}/classify_mags.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def VirMatch(self):
        '''
        Predict the hosts for viral contigs using VirMatcher.
        '''
        wkdir=f'{self.wkdir}/virmatcher'
        cmd=[self.envs]
        cmd.extend(
            ['VirMatcher --preparer','--gtdbtk-out-dir',self.tredir,
            '--gtdbtk-in-dir',self.host_mags,'-v',self.votus,
            '-o',wkdir,'--threads 32 --python-aggregator']
        )
        shell=f'{self.wkdir}/virmatcher.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
