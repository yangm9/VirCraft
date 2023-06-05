from os import path
from ..general import utils
from ..identify.viridsop import VirScan
from ..identify.multiFind import MultiTools

class VirRef(VirScan):
    '''
    Generate vOTUs
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
    def cluster(self):
        '''
        Cluster the sequence and remove redundancy for FastA file.
        '''
        votus=f'{self.outdir}/{self.name}_votus.fa'
        cmd=['cd-hit-est','-i',self.fasta,'-o',votus,'-T',self.threads,
            '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000\n']
        return cmd,votus
    def votuQC(self,votus,cutoff=1500):
        cmd=[utils.selectENV('viral-id-sop')]
        tmp_cmd,__=self.checkv(votus)
        cmd.extend(tmp_cmd)
        wkdir=f'{self.outdir}/stat'
        ckdir=f'{self.outdir}/checkv'
        checkv_qual=f'{ckdir}/checkv/quality_summary.tsv'
        cmd.extend(['pie_plot.R',checkv_qual,'checkv_quality',wkdir])
        cmd.extend(['pie_plot.R',checkv_qual,'provirus',wkdir])
        cmd.extend(['pie_plot.R',checkv_qual,'miuvig_quality',wkdir])
        tmp_cmd=self.statFA(cutoff)
        cmd.extend(self.statFA(cutoff))
        cmd.extend(tmp_cmd)
        MT=MultiTools(
            fasta=votus,
            outdir=self.outdir,
            threads=self.threads
        )
        tmp_cmd,vbdir=MT.vibrant()
        cmd.extend(tmp_cmd)
        votus_prefix=f'{self.name}_votus'
        vb_vir_type=f'{self.outdir}/VIBRANT_{votus_prefix}/VIBRANT_results_{votus_prefix}/VIBRANT_genome_quality_{votus_prefix}.tsv'
        vb_vir_type_m=f'{self.outdir}/vb_genome_type.tsv'
        vb_ckv_xls=f'{self.outdir}/vb_ckv_qual_type.xls'
        cmd.extend(
            ["sed '1s/scaffold/contig_id/'",vb_vir_type,vb_vir_type_m,'\n',
            'linkTab.py',checkv_qual,vb_vir_type_m,
            'left contig_id',vb_ckv_xls,'\n',
            'stacked_multi_bar_plot.R',vb_ckv_xls,wkdir,'\n']
        )
        return cmd
    def RmDup(self,cutoff=1500,unrun=False):
        cmd=[self.envs]
        tmp_cmd,votus=self.cluster()
        cmd.extend(tmp_cmd)
        cmd.extend(self.votuQC(votus,cutoff))
        shell=f'{self.outdir}/{self.name}_votu.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun:results=utils.execute(cmd)
        return results
