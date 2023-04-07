import os
from ..general import utils
from ..identify.multiFind import MultiTools

class vIdentify(MultiTools):
    '''
    Main Scripts
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def vFilter(self):
        cmd=[utils.selectENV('VirCraft')]
        cmd.extend(['merge_ctg_list.py',self.outdir,'\n'])
        return cmd
    def Identify(self,cutoff):
        #VirSorter2
        cmd=[utils.selectENV('viral-id-sop')]
        tmp_cmd,wkdir=self.virsorter(self.fasta,0)
        cmd.extend(tmp_cmd)
        vs2_partial_ctgs=f'{wkdir}/final-viral-score.tsv'
        vs2_out_ctgs_fa=f'{wkdir}/vs2-pass1/final-viral-combined.fa'
        vs2_partial_ctgs_fa=f'{self.outdir}/vs2_partial_viral_ctgs.fa'
        #vibrant
        shell=f'{self.outdir}/{self.name}_vs2_ctg.sh'
        utils.printSH(shell,cmd)
        cmd,wkdir=self.vibrant()
        shell=f'{self.outdir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell,cmd)
        vb_partial_ctgs_tab=f'{wkdir}/VIBRANT_results_{self.name}/VIBRANT_integrated_prophage_coordinates_{self.name}.tsv'
        vb_partial_ctgs_filt_tab=f'{wkdir}/VIBRANT_results_{self.name}/VIBRANT_integrated_prophage_coordinates_{self.name}.filt.tsv'
        partial_ctg_region=f'{wkdir}/partial_ctg_regions.bed'
        vb_partial_ctgs_fa=f'{self.outdir}/vb_partial_viral_ctgs.fa'
        #deepvirfinder
        cutoff=str(cutoff)
        cmd,wkdir=self.deepvirfinder(cutoff)
        shell=f'{self.outdir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell,cmd)
        
        #multiple run
        cmd=[utils.selectENV('')]
        cmd.extend(['multithreads.pl',self.outdir,'ctg.sh 3\n'])
        cmd.extend(self.vFilter())
        #full ctgs
        full_ctgs_li=f'{self.outdir}/full_viral_ctgs.list'
        full_ctgs_fa=f'{wkdir}/full_viral_ctgs.fa'
        cmd.extend(
            ['extrSeqByName.pl',full_ctgs_li,self.fasta,full_ctgs_fa,'\n']
        )
        #vs2 partial
        vs2_partial_ctgs_li=f'{self.outdir}/vs2_partial_viral_ctgs.list'
        cmd.extend(
            ['extrSeqByName.pl',vs2_partial_ctgs_li,
            vs2_out_ctgs_fa,vs2_partial_ctgs_fa,'\n']
        )
        #vb partial
        vb_partial_ctgs_li=f'{self.outdir}/vb_partial_viral_ctgs.list'
        cmd.extend(
            ["sed -i '1i\\fragment'",vb_partial_ctgs_li,'\n',
            'linkTab.py',vb_partial_ctgs_li,vb_partial_ctgs_tab,
            'left fragment',vb_partial_ctgs_filt_tab,'\n',
            'cut -f 1,6,7',vb_partial_ctgs_filt_tab,
            "|sed '1d' >",partial_ctg_region,'\n',
            'bedtools getfasta -fi',self.fasta,'-bed',partial_ctg_region,
            '-fo',vb_partial_ctgs_fa,'\n']
        )
        #cat all
        all_viral_ctgs=f'{self.outdir}/all_viral_ctgs.fa'
        cmd.extend(
            ['cat',full_ctgs_fa,vs2_partial_ctgs_fa,vb_partial_ctgs_fa,
            '>',all_viral_ctgs,'\n']
        )
        cmd.extend(tmp_cmd)
        shell=f'{self.outdir}/{self.name}_find_vir.sh'
        utils.printSH(shell,cmd)
        results=0#utils.execute(cmd)
        return results
