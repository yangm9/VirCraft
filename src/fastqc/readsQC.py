import os
from ..general import cmdExec,general
from ..config.cfgInfo import VirCfg

class Reads(VirCfg):
    '''
    FastQ processing class.
    '''
    envs=general.selectENV('VirCraft')
    def __init__(self,fq1='',fq2='',outdir='',*args,**kwargs):
        super().__init__()
        self.fastqs=[fq1,fq2]
        self.basename_fq1=os.path.basename(self.fastqs[0])
        self.basename_fq2=os.path.basename(self.fastqs[1])
        self.outdir=os.path.abspath(outdir)
        self.samp=self.basename_fq1.replace('_1.fastq','')
        self.samp=self.samp.replace('_1.fq','')
        general.mkdir(self.outdir)

class QualCtrl(Reads):
    def __init__(self,fq1='',fq2='',outdir='',process='',*args,**kwargs):
        super().__init__(fq1,fq2,outdir,*args,**kwargs)
        self.process=process
    def filtReads(self):
        wkdir=f'{self.outdir}/fastp'
        out_fq1=f'{wkdir}/{self.basename_fq1}'
        out_fq2=f'{wkdir}/{self.basename_fq2}'
        filt_rpts=f'{wkdir}/{self.samp}_report.html'
        out_fq_list=f'{wkdir}/{self.samp}_list.txt'
        general.mkdir(wkdir)
        cmd=''
        if 'f' in self.process:
            cmd=['fastp','-i',self.fastqs[0],'-o',out_fq1,
                '-I',self.fastqs[1],'-O',out_fq2,
                '-w 16 -q 20 -u 20 -g -c -W 5 -3 -l 50','-h',filt_rpts,'\n',
                'ls',out_fq1,out_fq2,'>',out_fq_list,'\n']
        else:
            cmd=['ls',self.fastqs[0],self.fastqs[1],'>',out_fq_list,'\n']
        return cmd,out_fq_list
    def fastUniq(self,fq_list):
        wkdir=f'{self.outdir}/fastuniq'
        print('======'+self.samp+'\n')
        in_fq_list=f'{self.outdir}/fastp/{self.samp}_list.txt'
        out_fq1=f'{wkdir}/{self.basename_fq1}'
        out_fq2=f'{wkdir}/{self.basename_fq2}'
        general.mkdir(wkdir)
        cmd=['fastuniq','-i',in_fq_list,'-o',out_fq1,'-p',out_fq2,'\n']
        fastqs=[out_fq1,out_fq2]
        return cmd,fastqs
    def decontaminate(self,fq1,fq2):
        wkdir=f'{self.outdir}/decontaminate'
        prefix=f'{wkdir}/{self.samp}'
        contam_db=self.confDict['ContamDB']
        cmd=['bowtie2','-p 40 -N 1','-x',contam_db,
            '-l',fq1,'-2',fq2,'--un-conc',prefix,'\n']
        fastqs=[f'{prefix}_1.fq',f'{prefix}_2.fq']
        return cmd,fastqs
    def readqc(self):
        cmd=[self.envs]
        tmp_cmd,fq_list=self.filtReads()
        cmd.extend(tmp_cmd)
        fastqs=[]
        if 'u' in self.process:
            tmp_cmd,fastqs=self.fastUniq(fq_list)
            cmd.extend(tmp_cmd)
        if 'c' in self.process:
            tmp_cmd,fastqs=self.decontaminate(fastqs[0],fastqs[1])
            cmd.extend(tmp_cmd)
        cmd.extend(
            [f'ln -s {fastqs[0]} {self.outdir}/{self.basename_fq1}\n',
            f'ln -s {fastqs[1]} {self.outdir}/{self.basename_fq2}\n']
        )
        shell=f'{self.outdir}/{self.samp}_readsqc.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
