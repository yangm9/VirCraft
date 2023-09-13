import os
from ..general import utils
from .bioseq import Reads

class QualCtrl(Reads):
    "FastQ cuality control class."
    envs=utils.selectENV('VC-ReadsQC')
    def __init__(self,fq1='',fq2='',outdir='',threads=8,*args,**kwargs):
        super().__init__(fq1,fq2,outdir,*args,**kwargs)
        self.threads=str(threads)
    def filtReads(self,process):
        wkdir=f'{self.outdir}/fastp'
        out_fq1=f'{wkdir}/{self.basename_fq1}'.replace('.gz','')
        out_fq2=f'{wkdir}/{self.basename_fq2}'.replace('.gz','')
        filt_rpts=f'{wkdir}/{self.samp}_report.html'
        out_fq_list=f'{wkdir}/{self.samp}_list.txt'
        utils.mkdir(wkdir)
        cmd=''
        out_fastqs=[]
        if 'f' in process:
            cmd=['fastp','-i',self.fastqs[0],'-o',out_fq1,
                '-I',self.fastqs[1],'-O',out_fq2,'-w',self.threads,
                self.confDict['FastpOpts'],'-h',filt_rpts,'\n',
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
        utils.mkdir(wkdir)
        cmd=['fastuniq','-i',in_fq_list,'-o',out_fq1,'-p',out_fq2,'\n']
        fastqs=[out_fq1,out_fq2]
        return cmd,fastqs
    def decontaminate(self,fq1,fq2):
        wkdir=f'{self.outdir}/decontaminate'
        prefix=f'{wkdir}/{self.samp}'
        contam_db=self.confDict['ContamDB']
        utils.mkdir(wkdir)
        outfqs=[f'{prefix}.1',f'{prefix}.2']
        fastqs=[f'{prefix}_1.fq',f'{prefix}_2.fq']
        cmd=['bowtie2','-p',self.threads,'-N 1','-x',contam_db,
            '-1',fq1,'-2',fq2,'--un-conc',prefix,'>/dev/null\n',
            'mv',outfqs[0],fastqs[0],'\n','mv',outfqs[1],fastqs[1],'\n']
        return cmd,fastqs
    def readqc(self,process='fuc',unrun=False,clear=False):
        cmd=[self.envs]
        fastqs=self.fastqs
        unused_fqs=[]
        if 'f' in process:
            tmp_cmd,fastqs,fq_list=self.filtReads(process)
            cmd.extend(tmp_cmd)
            unused_fqs.extend(fastqs)
        if 'u' in process:
            tmp_cmd,fastqs=self.fastUniq(fq_list)
            cmd.extend(tmp_cmd)
            unused_fqs.extend(fastqs)
        if 'c' in process:
            tmp_cmd,fastqs=self.decontaminate(fastqs[0],fastqs[1])
            cmd.extend(tmp_cmd)
            unused_fqs.extend(fastqs)
        lnk_fq1=os.path.basename(fastqs[0])
        lnk_fq2=os.path.basename(fastqs[1])
        cmd.extend(
            [f'ln -s {fastqs[0]} {self.outdir}/{lnk_fq1}\n',
            f'ln -s {fastqs[1]} {self.outdir}/{lnk_fq2}\n']
        )
        if clear:
            unused_fqs.remove(fastqs[0])
            unused_fqs.remove(fastqs[1])
            unused_fqs_str=' '.join(unused_fqs)
            cmd.extend(['rm -f',unused_fqs_str,'\n'])
        shell=f'{self.outdir}/{self.samp}_readsqc.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
