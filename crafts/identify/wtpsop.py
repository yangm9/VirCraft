import os
from ..general import cmdExec
from ..general import general
from ..config.config import Seq

class posiWTP(Seq):
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
        general.mkdir(wkdir)
        cmd=['nextflow run ~/.nextflow/assets/replikation/What_the_Phage',
            '--fasta',self.fasta,'--databases',self.confDict['wtpDB'],
            '--cachedir',self.confDict['wtpIMG'],'--output',resdir,
            '--workdir',tmpdir,'--cores',self.threads,
            '-profile local,singularity --filter 10000\n']
        return cmd,wkdir
    def checkv(self,in_fa:str):
        wkdir=f'{self.outdir}/checkv'
        general.mkdir(wkdir)
        cmd=['checkv','end_to_end',in_fa,wkdir,
             '-d',self.confDict['CheckvDB'],
             '-t',self.threads,'\n']
        provir_fna=f'{wkdir}/proviruses.fna'
        vir_fna=f'{wkdir}/viruses.fna'
        out_fa=f'{wkdir}/combined.fna'
        cmd.extend(
            ['cat',provir_fna,vir_fna,'>',out_fa,'\n']
        )
        return cmd,out_fa
    def annotate(self,indir:str):
        vs2_fa=f'{indir}/for-dramv/final-viral-combined-for-dramv.fa'
        vs2_tab=f'{indir}/for-dramv/viral-affi-contigs-for-dramv.tab'
        wkdir=f'{self.outdir}/dramv-annotate'
        cmd=['DRAM-v.py annotate','-i',vs2_fa,'-v',vs2_tab,'-o',wkdir,
            '--threads',self.threads,'--skip_trnascan --min_contig_size 1000\n']
        dramv_annot=f'{wkdir}/annotations.tsv'
        wkdir=f'{self.outdir}/dramv-distill'
        cmd.extend(['DRAM-v.py distill','-i',dramv_annot,'-o',wkdir,'\n'])
        return cmd
    def curate(self):
        checkv_dir=f'{self.outdir}/checkv'
        wkdir=f'{self.outdir}/curation'
        general.mkdir(wkdir)
        vir_score=f'{self.outdir}/vs2-pass1/final-viral-score.tsv'
        curation_score=f'{wkdir}/final-viral-score.tsv'
        contamination=f'{checkv_dir}/contamination.tsv'
        cmd=["sed '1s/seqname/contig_id/'",vir_score,'>',curation_score,'\n']
        cura_vs2_chkv=f'{wkdir}/curation_vs2_checkv.tsv'
        cmd.extend(
            ['linkTab.py',curation_score,contamination,
            'left contig_id',cura_vs2_chkv,'\n',
            'vCurator.py',self.outdir,'\n']
        )
        combined_fna=f'{checkv_dir}/combined.fna'
        combined_modi_fna=f'{checkv_dir}/combined_modi.fna'
        curated_contigs_xls=f'{wkdir}/curated_contigs.xls'
        curated_contigs_list=f'{wkdir}/curated_contigs_id.list'
        virus_posi_fna=f'{wkdir}/virus_positive.fna'
        cmd.extend(
            ['cut -f 1',curated_contigs_xls,
            '|grep -v "contig_id" >',curated_contigs_list,'\n',
            "sed 's/_[12] / /'",combined_fna,'>',combined_modi_fna,'\n',
            'extrSeqByName.pl',curated_contigs_list,combined_modi_fna,
            virus_posi_fna,'\n']
        )
        return cmd
    def Identify(self):
        cmd=[self.envs]
        #Step 1 Run VirSorter2
        tmp_cmd,wkdir=self.virsorter(self.fasta)
        cmd.extend(tmp_cmd)
        #Step 2 Run CheckV
        vs2_fa=f'{wkdir}/final-viral-combined.fa'
        tmp_cmd,checkv_fa=self.checkv(vs2_fa)
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
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
