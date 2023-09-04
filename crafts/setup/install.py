from ..general import utils

class deploy:
    '''
    Install the necessary softwares and conda environments for VirCraft.
    '''
    env_create='mamba create -c bioconda -c conda-forge'
    def __init__(self,outdir=''):
        self.outdir=outdir
        utils.mkdir(self.outdir)
    def readsqcENV(self):
        wkdir=f'{self.outdir}/VC-ReadsQC'
        cmd=[self.env_create,'-p',wkdir,'fastp=0.23.2 fastuniq=1.1 bowtie2=2.4.4 python=3.9.16 sortmerna=4.3.4 -y >/dev/null 2>&1']
        return cmd
    def assemblyENV(self):
        wkdir=f'{self.outdir}/VC-Assembly'
        cmd=[self.env_create,'-p',wkdir,'megahit=1.2.9 spades=3.15.4 bwa=0.7.17 -y >/dev/null 2>&1']
        return cmd
    def vbEnv(self):
        wkdir=f'{self.outdir}/VC-VIBRANT'
        utils.mkdir(vbdir)
        cmd=[self.env_create,'-p',wkdir,'python=3.7 vibrant=1.2.1 scikit-learn=0.21.3 biopython -y >/dev/null 2>&1']
        return cmd
    def vs2Env(self):
        wkdir=f'{self.outdir}/VC-VirSorter2'
        cmd=[self.env_create,'-p',wkdir,'python=3.8 virsorter=2.2.4 numpy=1.20.0 -y >/dev/null 2>&1']
        return cmd
    def dvfEnv(self):
        wkdir=f'{self.outdir}/VC-DeepVirFinder'
        cmd=[self.env_create,'-p',wkdir,'deepvirfinder=2020.11.21 -y >/dev/null 2>&1']
        return cmd
    def ckvENV(self):
        wkdir=f'{self.outdir}/VC-CheckV'
        cmd=[self.env_create,'-p',wkdir,'python=3.8 checkv=1.0.1 diamond=2.0.15 hmmer=3.3.2 prodigal -y >/dev/null 2>&1']
        return cmd
    def gtdbtkENV(self):
        wkdir=f'{self.outdir}/VC-GTDBTK'
        cmd=[self.env_create,'-p',wkdir,'gtdbtk=2.1.1 numpy=1.20.0 -y >/dev/null 2>&1']
        return cmd
    def vContact2ENV(self):
        wkdir=f'{self.outdir}/VC-vContact2'
        cmd=[self.env_create,'-p',wkdir,'python=3.7 vcontact2=0.11.0 pytables biopython networkx numpy=1.19.0 pandas=0.25.3 scipy=1.6.1 scikit-learn=0.24.1 psutil pyparsing hdf5 clusterone mcl blast diamond=2.0.15 -y >/dev/null 2>&1\n',
        'wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -q\n',
        'mv cluster_one-1.0.jar',wkdir,'\n']
        return cmd
    def ksENV(self):
        wkdir=f'{self.outdir}/VC-KofamScan'
        cmd=[self.env_create,'-c r','-p',wkdir,'kofamscan=1.3.0,hmmer=3.3.2 -y >/dev/null 2>&1\n'
        return cmd
        


