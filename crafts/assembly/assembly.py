from ..general import utils
from ..data.bioseq import Reads
from ..data.bioseq import Seq

class Assembly(Reads):
    '''
    Assembly the clean reads combined spades and megahit.
    '''
    envs=utils.selectENV('VC-Assembly')
    def __init__(self,fq1='',fq2='',outdir='',threads=8):
        super().__init__(fq1,fq2,outdir)
        self.threads=str(threads)
        self.methDict={
            'm':self.megahit,
            's':self.spades
        }
    def spades(self,fastqs:list):
        '''
        Assemble metagenome by SPAdes.
        '''
        wkdir=f'{self.wkdir}/spades'
        utils.mkdir(wkdir)
        if len(fastqs)==1:
            input_para=f'-s {fastqs[0]}'
        elif len(fastqs)==2:
            input_para=f'--pe1-1 {fastqs[0]} --pe1-2 {fastqs[1]}'
        elif len(fastqs)==3:
            input_para=f'--pe1-1 {fastqs[0]} --pe1-2 {fastqs[1]} -s {fastqs[2]}'
        else:
            pass
        cmd=['spades.py',input_para,'-t',self.threads,'-o',wkdir,self.confDict['SPAdesOpts'],'\n']
        scaf=f'{wkdir}/scaffolds.fasta'
        return cmd,scaf
    def megahit(self,fastqs:list):
        '''
        Assemble metagenome by megahit.
        '''
        input_para=''
        other_paras=''
        wkdir=f'{self.wkdir}/megahit'
        tmpdir=f'{self.wkdir}/megahit.tmp'
        utils.mkdir(tmpdir)
        if len(fastqs)==1:
            input_para=f'-r {fastqs[0]}'
        elif len(fastqs)==2:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]}'
            other_paras='--continue'
        elif len(fastqs)==3:
            input_para=f'-1 {fastqs[0]} -2 {fastqs[1]} -r {fastqs[2]}'
            other_paras='--continue'
        else:
            pass
        cmd=['megahit',input_para,'-o',wkdir,
            '-t',self.threads,self.confDict['MegahitOpts'],
            '--tmp-dir',tmpdir,other_paras,'\n']
        scaf=f'{wkdir}/final.contigs.fa'
        return cmd,scaf
    def unmapReads(self,scaf:str):
        '''
        Align the FastQs back to Assembled Contigs.
        '''
        wkdir=f'{self.wkdir}/alignment'
        utils.mkdir(wkdir)
        bwa_idx=f'{wkdir}/scaffoldsIDX'
        unused_bam=f'{wkdir}/unused_reads.bam'
        unused_fq_1=f'{wkdir}/unused_reads_1.fq'
        unused_fq_2=f'{wkdir}/unused_reads_2.fq'
        unused_fq_s=f'{wkdir}/unused_reads_s.fq'
        cmd=['bwa index -a bwtsw',scaf,'-p',bwa_idx,'\n',
            'bwa mem','-t',self.threads,bwa_idx,self.fastqs[0],self.fastqs[1],
            '|samtools view -b >',unused_bam,'\n',
            'samtools fastq -N',unused_bam,'-1',unused_fq_1,'-2',unused_fq_2,
            '-s',unused_fq_s,'\n']
        unused_fqs=[unused_fq_1,unused_fq_2,unused_fq_s]
        return cmd,unused_fqs
    def mixAsse(self,fastqs,process='m'):
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
            tmp_cmd,unused_fqs=self.unmapReads(tmp_scaf)
            cmd.extend(tmp_cmd)
            tmp_cmd,tmp_scaf=self.methDict[process[1]](unused_fqs)
            cmd.extend(tmp_cmd)
            scafs.append(tmp_scaf)
            cmd.extend(['cat',scafs[0],scafs[1],'>',final_scaf,'\n'])
        elif steps==1:
            cmd.extend(['ln -s',tmp_scaf,final_scaf,'\n'])
        else: 
            pass
        return cmd,final_scaf
    def Assemble(self,process='m',cutoff=1500,unrun=False,clear=False):
        cmd=[self.envs]
        scafs=[]
        tmp_cmd,scaf=self.mixAsse(self.fastqs,process)
        cmd.extend(tmp_cmd)
        FastA=Seq(scaf,self.outdir)
        cmd.extend(FastA.statFA(cutoff))
        if clear and len(process)==2:
            alndir=f'{self.wkdir}/alignment'
            scafIdx=f'{alndir}/scaffoldsIDX*'
            cmd.extend(['rm -f',scafIdx,'\n'])
        shell=f'{self.shelldir}/{self.samp}_assembly.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(shell)
        return results
