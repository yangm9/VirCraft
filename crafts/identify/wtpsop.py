import os
from ..general import utils
from ..data.bioseq import VirSeq
from ..data.bioseq import ORF

class posiWTP(VirSeq):
    '''
    posiWTP: positive WTP
    According to 10.1038/s42003-022-04027-y
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def wtp(self):
        wkdir=f'{self.outdir}/1.wtp'
        resdir=f'{wkdir}/result'
        tmpdir=f'{wkdir}/tempt'
        utils.mkdir(wkdir)
        cmd=['nextflow run ~/.nextflow/assets/replikation/What_the_Phage',
            '--fasta',self.fasta,'--databases',self.confDict['wtpDB'],
            '--cachedir',self.confDict['wtpIMG'],'--output',resdir,
            '--workdir',tmpdir,'--cores',self.threads,
            '-profile local,singularity --filter 10000\n']
        return cmd,resdir
    def curate(self):
        checkv_dir=f'{self.outdir}/checkv'
        wkdir=f'{self.outdir}/curation'
        wtpdir=f'{self.outdir}/wtp/{self.name}'
        checkvdir=f'{self.outdir}/checkv'
        posi_ctg_list=f'{wkdir}/{self.name}_posi_ctg.list'
        utils.mkdir(wkdir)
        orf_f='{self.outdir}/prodigal/{self.name}.orf.faa'
        eggout='{self.outdir}/eggnog/{self.name}.eggout.emapper.annotations'
        keywords_fa='{wkdir}/{self.name}_vkeywd_posi_ctg.fa'
        posi_ctg_fa='{wkdir}/{self.name}_v_posi_ctg.fa'
        cmd.extend(
            ['virus_identity.py',self.fasta,orf_f,keywords_fa,'\n',
            'filtByViralPredictionTools.py',wtpdir,checkvdir,
            '>',posi_ctg_list,'\n',
            'extrSeqByName.pl',posi_ctg_list,keywords_fa,posi_ctg_fa,'\n']
        )
        return cmd
    def Identify(self):
        cmd=[self.envs]
        #Step 1 Run WTP
        tmp_cmd,resdir=self.wtp(self.fasta)
        cmd.extend(tmp_cmd)
        #Step 2 Predict Genes
        self.fasta=f'{resdir}/phage_positive_contigs/{self.name}_positive_contigs.fa'
        tmp_cmd,orf_faa=self.genePred()
        cmd.extend(tmp_cmd)
        #Step 3 eggnog annotation
        SeqO=ORF(orf=orf_faa,outdir=self.outdir)
        cmd.extend(SeqO.eggnogAnno())
        tmp_cmd,checkv_fa=self.checkv(wtp_fa)
        cmd.extend(tmp_cmd)
        #Step 3 rerun VirSorter2 
        tmp_cmd,wkdir=self.virsorter(checkv_fa)
        cmd.extend(tmp_cmd)
        #Step 4 DRAMv annotation
        tmp_cmd=self.annotate(wkdir) 
        cmd.extend(tmp_cmd)
        #Step 5 Curation
        tmp_cmd=self.curate()
        cmd.extend(tmp_cmd)
        #Generate shell and exeute it
        shell=f'{self.outdir}/{self.name}_find_vir.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(shell)
        return results
