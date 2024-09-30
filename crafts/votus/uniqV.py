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
    def cdhit_cluster(self):
        '''
        Cluster the sequence and remove redundancy for FastA file using CD-HIT.
        '''
        votus=f'{self.outdir}/{self.name}_votus.fa'
        cmd=[utils.selectENV('VC-General')]
        cmd=['cd-hit-est','-i',self.fasta,'-o',votus,'-T',self.threads,
            '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000\n']
        return cmd,votus
    def blast_cluster(self):
        '''
        Rapid genome clustering based on pairwise ANI provided by CheckV.
        '''
        blastdb=f'{self.outdir}/{self.name}.db'
        blast_out=f'{self.outdir}/{self.name}.blast'
        ani_out=f'{self.outdir}/{self.name}.ani'
        clusters=f'{self.outdir}/{self.name}.clusters'
        votu_list=f'{self.outdir}/{self.name}.votulist'
        votus=f'{self.outdir}/{self.name}_votus.fa'
        cmd=['makeblastdb','-in',self.fasta,'-dbtype nucl','-out',blastdb,'\n',
            'blastn','-query',self.fasta,'-db',blastdb,'-num_threads',self.threads,
            '-out',blast_out,"-outfmt '6 std qlen slen' -max_target_seqs 10000\n",
            'anicalc.py','-i',blast_out,'-o',ani_out,'\n',
            'aniclust.py','--fna',self.fasta,'--ani',ani_out,'--out',clusters,
            '--min_ani 95 --min_tcov 85 --min_qcov 0\n',
            'cut -f 1',clusters,'>',votu_list,'\n',
            'extrSeqByName.pl',votu_list,self.fasta,votus,'\n']
        return cmd,votus
    def cluster(self,method):
        if method=='blast':
            return self.blast_cluster()
        elif method=='cdhit':
            return self.cdhit_cluster()
        else:
            raise Exception('method should only be "blast" or "cdhit".')
    def votuQC(self,votus,cutoff=1500):
        cmd,__=self.checkv(votus)
        wkdir=f'{self.outdir}/stat'
        ckdir=f'{self.outdir}/checkv'
        checkv_qual=f'{ckdir}/quality_summary.tsv'
        cmd.append(utils.selectENV('VC-General'))
        cmd.extend(['pie_plot.R',checkv_qual,'checkv_quality',wkdir,'\n'])
        cmd.extend(['pie_plot.R',checkv_qual,'provirus',wkdir,'\n'])
        cmd.extend(['pie_plot.R',checkv_qual,'miuvig_quality',wkdir,'\n'])
        tmp_cmd=self.statFA(cutoff)
        cmd.extend(self.statFA(cutoff))
        cmd.extend(tmp_cmd)
        MT=MultiTools(
            fasta=votus,
            outdir=self.outdir,
            threads=self.threads
        )
        tmp_cmd,vbdir=MT.vibrant(str(cutoff))
        cmd.extend(tmp_cmd)
        votus_prefix=f'{self.name}_votus'
        vb_vir_info=f'{self.outdir}/VIBRANT_{votus_prefix}/VIBRANT_results_{votus_prefix}/VIBRANT_genome_quality_{votus_prefix}.tsv'
        vb_ckv_tsv=f'{wkdir}/votus_lifetype_quality.tsv'
        cmd.append(utils.selectENV('VC-General'))
        cmd.extend(
            ['votus_lifetype_quality.py',checkv_qual,vb_vir_info,vb_ckv_tsv,'\n',
            'votus_lifetype_quality_barplot.R',vb_ckv_tsv,wkdir,'\n']
        )
        return cmd
    def binning():
        cmd=[utils.selectENV('VC-General')]
        wkdir=f'{self.outdir}'
        
        return cmd
    def RmDup(self,cutoff=1500,unrun=False,method='blast'):
        cmd=[self.envs]
        tmp_cmd,votus=self.cluster(method)
        cmd.extend(tmp_cmd)
        cmd.extend(self.votuQC(votus,cutoff))
        shell=f'{self.outdir}/{self.name}_votu.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun:results=utils.execute(shell)
        return results
