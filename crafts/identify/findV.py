import os
from ..general import utils
from ..identify.multiFind import MultiTools

class vIdentify(MultiTools):
    '''
    Main Scripts
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(int(threads)//self.BATCH_SIZE*2)
    def vFilter(self):
        #merge Contigs
        cmd=[utils.selectENV('VirCraft')]
        full_ctgs_li=f'{self.outdir}/full_viral_ctgs.list'
        full_ctgs_fa=f'{self.outdir}/full_viral_ctgs.fa'
        score_xls=f'{self.outdir}/all_viral_ctgs.score.xls'
        score_filt_xls=utils.insLable(score_xls,'gt2')
        score_filt_ctgs_li=f'{self.outdir}/viral_ctgs_filt.list'
        cmd.extend(
            ['merge_ctg_list.py',self.name,self.outdir,'\n',
            "awk -F '\\t' 'NR==1 || $20>=2'",score_xls,'>',score_filt_xls,'\n',
            'cut -f 1',score_filt_xls,"|sed '1d' >",score_filt_ctgs_li,'\n',
            'extrSeqByName.pl',full_ctgs_li,self.fasta,full_ctgs_fa,'\n']
        )
        #vs2 partial
        wkdir=f'{self.outdir}/vs2-pass1'
        vs2_partial_ctgs_li=f'{self.outdir}/vs2_partial_viral_ctgs.list'
        vs2_out_ctgs_fa=f'{wkdir}/final-viral-combined.fa'
        vs2_partial_ctgs_fa=f'{self.outdir}/vs2_partial_viral_ctgs.fa'
        cmd.extend(
            ['extrSeqByName.pl',vs2_partial_ctgs_li,
            vs2_out_ctgs_fa,vs2_partial_ctgs_fa,'\n']
        )
        #vb partial
        wkdir=f'{self.outdir}/VIBRANT_{self.name}'
        vb_partial_ctgs_li=f'{self.outdir}/vb_partial_viral_ctgs.list'
        vb_partial_ctgs_tab=f'{wkdir}/VIBRANT_results_{self.name}/VIBRANT_integrated_prophage_coordinates_{self.name}.tsv'
        vb_partial_ctgs_filt_tab=f'{wkdir}/VIBRANT_results_{self.name}/VIBRANT_integrated_prophage_coordinates_{self.name}.filt.tsv'
        partial_ctg_region=f'{wkdir}/partial_ctg_regions.bed'
        vb_partial_ctgs_fa=f'{self.outdir}/vb_partial_viral_ctgs.fa'
        ins_head_cmd=f'''if [ -s "{vb_partial_ctgs_li}" ]
then
    sed -i '1i\\fragment' {vb_partial_ctgs_li}
else
    echo 'fragment' > {vb_partial_ctgs_li}
fi
'''
        cmd.extend(
            [ins_head_cmd,'linkTab.py',vb_partial_ctgs_li,vb_partial_ctgs_tab,
            'left fragment',vb_partial_ctgs_filt_tab,'\n',
            'cut -f 1,6,7',vb_partial_ctgs_filt_tab,
            "|sed '1d' >",partial_ctg_region,'\n',
            'bedtools getfasta -fi',self.fasta,'-bed',partial_ctg_region,
            '-fo',vb_partial_ctgs_fa,'\n']
        )
        #cat all
        all_viral_ctgs=f'{self.outdir}/all_viral_ctgs.fa'
        tool_filt_ctgs=f'{self.outdir}/all_viral_ctgs_score_gt2.fa'
        cmd.extend(
            ['cat',full_ctgs_fa,vs2_partial_ctgs_fa,
            vb_partial_ctgs_fa,'>',all_viral_ctgs,'\n',
            'extrSeqByName.pl',score_filt_ctgs_li,all_viral_ctgs,
            tool_filt_ctgs,'\n']
        )
        tmp_cmd,checkv_fa=self.checkv(tool_filt_ctgs)
        cmd.extend(tmp_cmd)
        return cmd
    def Identify(self,cutoff=1500,unrun=False):
        #VirSorter2
        cmd=[utils.selectENV('viral-id-sop')]
        tmp_cmd,wkdir=self.virsorter(self.fasta,0,cutoff)
        cmd.extend(tmp_cmd)
        #vibrant
        self.threads=str(int(self.threads)//2)
        shell=f'{self.outdir}/{self.name}_vs2_ctg.sh'
        utils.printSH(shell,cmd)
        cmd,wkdir=self.vibrant()
        shell=f'{self.outdir}/{self.name}_vb_ctg.sh'
        utils.printSH(shell,cmd)
        #deepvirfinder
        cutoff=str(cutoff)
        cmd,wkdir=self.deepvirfinder(cutoff)
        shell=f'{self.outdir}/{self.name}_dvf_ctg.sh'
        utils.printSH(shell,cmd)
        #multiple run
        cmd=[utils.selectENV('')]
        cmd.extend(['multithreads.pl',self.outdir,'ctg.sh 3\n'])
        cmd.extend(self.vFilter())
        shell=f'{self.outdir}/{self.name}_find_vir.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(cmd)
        return results
