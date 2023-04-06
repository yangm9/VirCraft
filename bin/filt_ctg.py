#!/usr/bin/env python3
import pandas as pd
from shutil import which

PathDict={
    'virsorter2'='vs2-pass1/final-viral-score.tsv',
    'vibrant':'VIBRANT_scaffolds/VIBRANT_results_scaffolds/VIBRANT_machine_scaffolds.tsv',
    'deepvirfinder'='deepvirfinder'
}

FiltDict={
    'virsorter2':'dsDNAphage>=0.9',
    'vibrant':"prediction=='virus'",
    'deepvirfinder':'score > 0.9 and pvalue<=0.1'
}

NameDict={
    'virsorter2':'seqname',
    'vibrant':'scaffold',
    'deepvirfinder':'name'
}

FastaDict={
    'virsorter2':'vs2-pass1/final-viral-combined.fa',
    'vibrant':'',
    'deepvirfinder':''
}

def vCtgFilt(tool,wkdir):
    result=wkdir+'/'+PathDict[tool]
    if tool=='deepvirfinder':
        file_name=os.listdir(result)[0]
        result+='/'+file_name
    df=pd.read_csv(result,sep='\t')
    df=df.query(FiltDict[tool])
    return df

def vCtgMerge(wkdir):
    full_ctgs=[]
    partial_ctgs=[]
    for tool in FiltDict.keys():
        df=vCtgFilt(tool,wkdir)
        if tool=='virsorter2':
            df['seqname']=df['seqname'].apply(lambda x:x.rsplit('_',1)[0] if x.endswith('full') or x.endswith('lt2gene') else x)
            full_tmps=df.query('not seqname.str.contains("partial")')['seqname'].tolist()
            partial_tmps=df.query('seqname.str.contains("partial")')['seqname'].tolist()
            full_ctgs.extend(full_tmps)
            vs2_partial_ctgs.extend(partial_tmps)
        elif tool=='vibrant':
            full_tmps=df.query('not scaffold.str.contains("fragment")')['scaffold']
            partial_tmps=df.query('scaffold.str.contains("fragment")')['scaffold']
            full_ctgs.extend(full_tmps)
            vb_partial_ctgs.extend(partial_tmps)
        else:
            full_ctgs.extend(df[NameDict[tool]].tolist())
    full_ctgs=list(set(full_ctgs))
    partial_ctgs=list(set(partial_ctgs))
    return full_ctgs,vs2_partial_ctgs,vb_partial_ctgs

def listToFile(list_l,list_f):
    LIST=open(list_f,'w')
    for ctg in list_l:
        LIST.write(f'{ctg}\n')
    return 0

def selectENV(env:str):
    bin_dir=os.path.abspath(__file__)
    conda_path=which('conda')
    condash_path='/'.join(conda_path.split('/')[0:-2])
    condash_path+='/etc/profile.d/conda.sh'
    envs=''
    if os.path.exists(condash_path):
        envs=f'. "{condash_path}"\nconda activate {env}\n'
    else:
        conda_bin=os.path.dirname('conda_path')
        envs=f'export PATH="{conda_bin}:$PATH"\n'
    envs+=f'export PATH="{bin_dir}:$PATH"\n'
    return envs

def getVCtg(fasta,wkdir):
    full_ctgs,vs2_partial_ctgs,vb_partial_ctgs=vCtgMerge(wkdir)
    full_ctgs_li=f'{wkdir}/full_viral_ctgs.list'
    vs2_partial_ctgs_li=f'{wkdir}/vs2_partial_viral_ctgs.list'
    vb_partial_ctgs_li=f'{wkdir}/vb_partial_viral_ctgs.list'
    vb_partial_ctgs_tab=f'{wkdir}/VIBRANT_scaffolds/VIBRANT_results_scaffolds/VIBRANT_integrated_prophage_coordinates_scaffolds.tsv'
    vb_partial_ctgs_tab=f'{wkdir}/VIBRANT_scaffolds/VIBRANT_results_scaffolds/VIBRANT_integrated_prophage_coordinates_scaffolds.tsv'
    vb_partial_ctgs_filt_tab=f'{wkdir}/VIBRANT_scaffolds/VIBRANT_results_scaffolds/VIBRANT_integrated_prophage_coordinates_scaffolds.filt.tsv'
    full_ctgs_fa=f'{wkdir}/full_viral_ctgs.fa'
    vs2_out_ctgs_fa=f'{wkdir}/vs2-pass1/final-viral-combined.fa'
    vs2_partial_ctgs_fa=f'{wkdir}/vs2_partial_viral_ctgs.fa'
    vb_partial_ctgs_fa=f'{wkdir}/vb_partial_viral_ctgs.fa'
    listToFile(full_ctgs,full_ctgs_li)
    listToFile(vs2_partial_ctgs,vs2_partial_ctgs_li)
    listToFile(vb_partial_ctgs,vb_partial_ctgs_li)
    script_dir=os.path.abspath(__file__)
    cmd_txt=selectENVs('VirCraft')
    cmd_txt+=f'''
extrSeqByName.pl {full_ctgs_li} {fasta} {full_ctgs_fa}
extrSeqByName.pl {vs2_partial_ctgs_li} {vs2_out_ctgs_fa} {vs2_partial_ctgs_fa}
sed -i '1i\\fragment' {vb_partial_ctgs_li}
linkTab.py {vb_partial_ctgs_li} {vb_partial_ctgs_tab} left fragment {vb_partial_ctgs_filt_tab}
cut -f 1,6,7 vb_partial_ctgs_filt_tab|sed '1d' > {wkdir}/regions.bed
bedtools getfasta -fi {fasta} -bed {wkdir}/regions.bed -fo {vb_partial_ctgs_fa}
cat {full_ctgs_fa} {vs2_partial_ctgs_fa} {vb_partial_ctgs_fa} >all_viral_ctgs.fa
'''
    result=os.system(cmd_txt)
    return results

