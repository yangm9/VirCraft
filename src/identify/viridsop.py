import os
from ..general import cmdExec,general
from ..config.cfgInfo import VirCfg

class VirScan(VirCfg):
    '''
    According to the Viral sequence identification SOP with VirSorter2 (https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3)
    '''
    n=0
    vs2_subcmds=['--keep-original-seq','--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv']
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__()
        self.fasta=os.path.abspath(fasta)
        self.outdir=os.path.abspath(outdir)
        self.threads=str(threads)
        general.mkdir(self.outdir)
    def virsorter(self,in_fa:str):
        idx=str(self.n+1)
        wkdir=f'{self.outdir}/vs2-pass{idx}'
        general.mkdir(wkdir)
        cmd=['virsorter run',vs2_subcmds[self.n],'-i',in_fa,
            '-d',self.confDict['virsorter2DB'],'-w',wkdir,
            '--include-groups dsDNAphage,ssDNA','-j',self.threads,
            '--min-length 5000 --min-score 0.5 all\n']
        self.n+=1
        return cmd,wkdir
    def checkv(self,in_fa:str)
        wkdir=f'{self.outdir}/checkv'
        general.mkdir(wkdir)
        cmd=['checkv','end_to_end',in_fa,wkdir,
             '-d',self.confDict['CheckVDB'],'-t',self.threads]
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
        general.mkdir(wkdir)
        cmd=['DRAM-v.py annotate','-i',vs2_fa,'-v',vs2_tab,'-o',wkdir,
            '--threads',self.threads,'--skip_trnascan --min_contig_size 1000\n']
        dramv_annot=f'{wkdir}/annotations.tsv'
        wkdir=f'{self.outdir}/dramv-distill'
        general.mkdir(wkdir)
        cmd.extend(['DRAM-v.py distill','-i',dramv_annot,'-o',wkdir,'\n'])
        return cmd
    def curate(self):
        wkdir=f'{self.outdir}/curation'
        general.mkdir(wkdir)
        vir_score=f'{self.outdir}/vs2-pass1/final-viral-score.tsv'
        curation_score=f'{wkdir}/final-viral-score.tsv'
        contamination=f'{self.outdir}/checkv/contamination.tsv'
        cmd=["sed '1s/seqname/contig_id/'",vir_score,'>',curation_score,'\n']
        cura_vs2_chkv=f'{wkdir}/curation_vs2_checkv.tsv'
        cmd.extend(
            ['linkTab.py',curation_score,contamination,
            'left contig_id',cura_vs2_chkv,'\n',
            'vCurator.py',self.outdir,'\n']
        )
        combined_modi_fna=f'{self.outdir}/checkv/combined_modi.fna'
        contig_id_list=f'{wkdir}/contigs_id.list'
        virus_posi_fna=f'{wkdir}/virus_positive.fna'
        cmd.extend(
            ["sed 's/_1 / /'",combined_fna,'>',combined_modi_fna,'\n',
            'extrSeqByName.pl',contig_id_list,combined_modi_fna,
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
        tmp_cmd=self.curate(wkdir)
        cmd.extend(tmp_cmd)
        #Generate shell and exeute it
        shell=f'{self.outdir}/find_vir.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
