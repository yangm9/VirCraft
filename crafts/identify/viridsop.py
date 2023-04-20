import os
from ..general import utils
from ..data.bioseq import Seq

class VirScan(Seq):
    '''
    According to the Viral sequence identification SOP with VirSorter2 (https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3)
    '''
    vs2_subcmds=['--keep-original-seq','--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv']
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def virsorter(self,in_fa:str,n:int):
        idx=str(n+1)
        wkdir=f'{self.outdir}/vs2-pass{idx}'
        utils.mkdir(wkdir)
        cmd=['virsorter run',self.vs2_subcmds[n],'-i',in_fa,
            '-d',self.confDict['Virsorter2DB'],'-w',wkdir,
            '--include-groups dsDNAphage,ssDNA','-j',self.threads,
            '--min-length 5000 --min-score 0.5 all\n']
        return cmd,wkdir
    def checkv(self,in_fa:str):
        wkdir=f'{self.outdir}/checkv'
        utils.mkdir(wkdir)
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
        utils.mkdir(wkdir)
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
    def Identify(self,unrun=False):
        cmd=[utils.selectENV('viral-id-sop')]
        #Step 1 Run VirSorter2
        tmp_cmd,wkdir=self.virsorter(self.fasta,0)
        cmd.extend(tmp_cmd)
        #Step 2 Run CheckV
        vs2_fa=f'{wkdir}/final-viral-combined.fa'
        tmp_cmd,checkv_fa=self.checkv(vs2_fa)
        cmd.extend(tmp_cmd)
        #Step 3 rerun VirSorter2 
        tmp_cmd,wkdir=self.virsorter(checkv_fa,1)
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
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
