import os
from ..general import cmdExec,general
from ..config.config import Reads

class QualCtrl(Reads):
    "FastQ cuality control class."
    envs=general.selectENV('VirCraft')
    def __init__(self,fq1='',fq2='',outdir='',*args,**kwargs):
        super().__init__(fq1,fq2,outdir,*args,**kwargs)
    def filtReads(self,process):
        wkdir=f'{self.outdir}/fastp'
        out_fq1=f'{wkdir}/{self.basename_fq1}'.replace('.gz','')
        out_fq2=f'{wkdir}/{self.basename_fq2}'.replace('.gz','')
        filt_rpts=f'{wkdir}/{self.samp}_report.html'
        out_fq_list=f'{wkdir}/{self.samp}_list.txt'
        general.mkdir(wkdir)
        cmd=''
        out_fastqs=[]
        if 'f' in process:
            cmd=['fastp','-i',self.fastqs[0],'-o',out_fq1,
                '-I',self.fastqs[1],'-O',out_fq2,
                '-w 16 -q 20 -u 20 -g -c -W 5 -3 -l 50','-h',filt_rpts,'\n',
                'ls',out_fq1,out_fq2,'>',out_fq_list,'\n']
            out_fastqs=[out_fq1,out_fq2]
        else:
            cmd=['ls',self.fastqs[0],self.fastqs[1],'>',out_fq_list,'\n']
            out_fastqs=self.fastqs
        return cmd,out_fastqs,out_fq_list
    def fastUniq(self,fq_list):
        wkdir=f'{self.outdir}/fastuniq'
        in_fq_list=f'{self.outdir}/fastp/{self.samp}_list.txt'
        out_fq1=f'{wkdir}/{self.basename_fq1}'.replace('.gz','')
        out_fq2=f'{wkdir}/{self.basename_fq2}'.replace('.gz','')
        general.mkdir(wkdir)
        cmd=['fastuniq','-i',in_fq_list,'-o',out_fq1,'-p',out_fq2,'\n']
        fastqs=[out_fq1,out_fq2]
        return cmd,fastqs
    def decontaminate(self,fq1,fq2):
        wkdir=f'{self.outdir}/decontaminate'
        prefix=f'{wkdir}/{self.samp}'
        contam_db=self.confDict['ContamDB']
        general.mkdir(wkdir)
        cmd=['bowtie2','-p 40 -N 1','-x',contam_db,
            '-l',fq1,'-2',fq2,'--un-conc',prefix,'\n']
        fastqs=[f'{prefix}_1.fq',f'{prefix}_2.fq']
        return cmd,fastqs
    def readqc(self,process='fuc'):
        cmd=[self.envs]
        fastqs=[]
        tmp_cmd,fastqs,fq_list=self.filtReads(process)
        cmd.extend(tmp_cmd)
        if 'u' in process:
            tmp_cmd,fastqs=self.fastUniq(fq_list)
            cmd.extend(tmp_cmd)
        if 'c' in process:
            tmp_cmd,fastqs=self.decontaminate(fastqs[0],fastqs[1])
            cmd.extend(tmp_cmd)
        lnk_fq1=os.path.basename(fastqs[0])
        lnk_fq2=os.path.basename(fastqs[1])
        cmd.extend(
            [f'ln -s {fastqs[0]} {self.outdir}/{lnk_fq1}\n',
            f'ln -s {fastqs[1]} {self.outdir}/{lnk_fq2}\n']
        )
        shell=f'{self.outdir}/{self.samp}_readsqc.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
