from ..general import utils
from ..data.bioseq import Reads
from ..data.bioseq import Seq

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.envs=utils.selectENV('assembly')
        self.threads=str(threads)
        self.methDict={
            's':self.spades,
            'm':self.megahit
        }
    def spades(self,fastqs:list):
        '''
        Assemble metagenome by SPAdes.
        '''
        wkdir=f'{self.outdir}/spades'
        utils.mkdir(wkdir)
        cmd=[]
        if len(fastqs)==1:
            cmd=['spades.py','-s',fastqs[0]]
        elif len(fastqs)==2:
            cmd=['spades.py','--pe1-1',fastqs[0],'--pe1-2',fastqs[1]]
        else:
            pass
        cmd.extend(
            ['-t',self.threads,'-o',wkdir,
            '--careful -m 1300 -k 21,33,55,77,99,127','\n']
        )
        scaf=f'{wkdir}/scaffolds.fasta'
        return cmd,scaf
    def megahit(self,fastqs:list):
        '''
        Assemble metagenome by megahit.
        '''
        input_para=''
        other_paras=''
        wkdir=f'{self.outdir}/megahit'
        tmpdir=f'{self.outdir}/megahit.tmp'
        utils.mkdir(tmpdir)
        if len(fastqs)==1:
            input_para=f'-r {fastqs[0]}'
        elif len(fastqs)==2:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]}'
            other_paras='--continue'
        else:
            pass
        cmd=['megahit',input_para,'-o',wkdir,
            '-t',self.threads,'-m','80000000000',
            '--tmp-dir',tmpdir,other_paras,'\n']
        scaf=f'{wkdir}/final.contigs.fa'
        return cmd,scaf
    def unmapReads(self,scaf:str):
        '''
        Align the FastQs back to Assembled Contigs.
        '''
        wkdir=f'{self.outdir}/alignment'
        utils.mkdir(wkdir)
        bwa_idx=f'{wkdir}/scaffoldsIDX'
        unused_sam=f'{wkdir}/unused_reads.sam'
        unused_fq=f'{wkdir}/unused_reads.fq'
        cmd=['bwa index -a bwtsw',scaf,'-p',bwa_idx,'\n',
            'bwa mem','-t',self.threads,bwa_idx,self.fastqs[0],self.fastqs[1],
            '|grep -v NM:i:>',unused_sam,'\n',
            'sam_to_fastq.py',unused_sam,'>',unused_fq,'\n']
        return cmd,unused_fq
    def mixAsse(self,fastqs,process='sm'):
        '''
        Analysis the Assembly process accordding to the process set.
        '''
        cmd=[]
        scafs=[]
        tmp_cmd,tmp_scaf=self.methDict[process[0]](fastqs)
        cmd.extend(tmp_cmd)
        scafs.append(tmp_scaf)
        final_scaf=f'{self.outdir}/final_assembly.fasta'
        steps=len(process)
        if steps==2:
            tmp_cmd,unused_fq=self.unmapReads(tmp_scaf)
            cmd.extend(tmp_cmd)
            tmp_cmd,tmp_scaf=self.methDict[process[1]]([unused_fq])
            cmd.extend(tmp_cmd)
            scafs.append(tmp_scaf)
            cmd.extend(['cat',scafs[0],scafs[1],'>',final_scaf,'\n'])
        elif steps==1:
            cmd.extend(['ln -s',tmp_scaf,final_scaf,'\n'])
        return cmd,final_scaf
    def Assemble(self,process='sm',cutoff=5000,unrun=False):
        cmd=[self.envs]
        scafs=[]
        tmp_cmd,scaf=self.mixAsse(self.fastqs,process)
        cmd.extend(tmp_cmd)
        FastA=Seq(scaf,self.outdir)
        cmd.extend(FastA.statFA(cutoff))
        shell=f'{self.outdir}/{self.samp}_assembly.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
