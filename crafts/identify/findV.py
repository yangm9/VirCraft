import os
import sys
from ..general import utils
from ..identify.multiFind import MultiTools

class vIdentify(MultiTools):
    '''
    Main Scripts of identify module.
    self.BATCH_SIZE initialized as 4 in VirCfg class will be used to divide the input threads into 2*self.BATCH_SIZE portions, with 2 allocated to VirSorter2, and the other 2 allocated to VIBRANT and DeepVirFinder respectively.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.allthreads=threads
        self.threads=int(threads)//(self.BATCH_SIZE*2)
    def vFilter(self):
        #merge Contigs
        score_tsv=f'{self.outdir}/all_viral_ctgs.score.tsv'
        score_filt_tsv=utils.insLable(score_tsv,'gt2')
        filt_ctgs_list=f'{self.outdir}/viral_ctgs_filt.list'
        filt_viral_ctgs=f'{self.outdir}/filted_viral_ctgs.fa'
        cmd=[utils.selectENV('VC-General')]
        tmp_cmd,cat_dir=self.contig_annotation_tool(filt_viral_ctgs)
        cmd.extend(
            ['merge_ctg_list.py',self.name,self.outdir,'\n',
            "awk -F '\\t' 'NR==1 || $21>=2'",score_tsv,'>',score_filt_tsv,'\n',
            'cut -f 1',score_filt_tsv,"|sed '1d' >",filt_ctgs_li,'\n',
            'extrSeqByName.pl',filt_ctgs_li,self.fasta,filt_viral_ctgs,'\n']
        )
        tmp_cmd,checkv_fa=self.checkv(filt_viral_ctgs)
        cmd.extend(tmp_cmd)
        checkv_dir=os.path.dirname(checkv_fa)
        filt_fna_files=f'{checkv_dir}/*.filt.fna'
        filt_checkv_fa=f'{self.outdir}/viral_positive_ctgs.fna'
        cmd.extend([utils.selectENV('VC-General')])
        cmd.extend(
            ['vir_qual_filt.py',checkv_dir,'\n',
            'cat',filt_fna_files,'>',filt_checkv_fa,'\n']
        )
        return cmd
    def Identify(self,cutoff=1500,unrun=False):
        try:
            if int(self.allthreads)<8:
                raise ValueError('The threads number must not be less than 8!!!')
        except ValueError as e:
            print(f'ERROR: {e}')
            exit(1)
        results=''
        cutoff=str(cutoff)
        self.threads=str(self.threads)
        #vibrant
        cmd,wkdir=self.vibrant(cutoff)
        shell=f'{self.outdir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell,cmd)
        #deepvirfinder
        cmd,wkdir=self.deepvirfinder(cutoff)
        shell=f'{self.outdir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell,cmd)
        #VirSorter2
        self.threads=str(int(self.allthreads)-(int(self.threads)*2))
        cmd,wkdir=self.virsorter(self.fasta,0,cutoff)
        shell=f'{self.outdir}/{self.name}_vs2_ctg.sh'
        utils.printSH(shell,cmd)
        #multiple run
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(['multithreads.pl',self.outdir,'ctg.sh 3\n'])
        shell=f'{self.outdir}/{self.name}_multi_find.sh'
        utils.printSH(shell,cmd)
        if not unrun: results=utils.execute(shell) 
        self.threads=str(int(self.allthreads))
        cmd=self.vFilter()
        shell=f'{self.outdir}/{self.name}_merge_vir.sh'
        utils.printSH(shell,cmd)
        if not unrun: results+=utils.execute(shell)
        return results
