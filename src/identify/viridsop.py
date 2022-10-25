
from ..config import setVari,conf
from ..process import cmdExec

groups, confDict, sampDict = conf.prepInfo(config)

def findVir(grp: str, outdir: str):
    envs = setVari.selectENV('VirCraft')
    find_vir_cmd = [envs]
    identify_dir = f'{outdir}/02.identify/{grp}'
    if not os.path.exists(identify_dir): os.makedirs(identify_dir)
    vs2_pass1_dir = f'{identify_dir}/vs2-pass1'
    if not os.path.exists(vs2_pass1_dir): os.makedirs(vs2_pass1_dir)
    
    scaffolds = f'{outdir}/01.assembly/{grp}/scaffolds.filt.gt5000.fasta'
    find_vir_cmd.extend(
        ['virsorter', 'run', '--keep-original-seq', '-i', scaffolds, 
         '-d ', confDict['virsorter2DB'], '-v', vs2_pass1_dir, 
         '--include-groups dsDNAphage,ssDNA', 
         '--min-length 0 --min-score 0.5 -j 28 all\n']
    )

    vs2_pass1_fasta=f'{vs2_pass0}/final-viral-combined.fa'
    checkv_dir=f'{identify_dir}/checkv'
    if not os.path.exists(checkv_dir): os.makedirs(checkv_dir)
    find_vir_cmd.extend(
        ['checkv', 'end_to_end', vs2_pass1_fasta, checkv_dir,
         '-d', confDict['CheckVDB'], '-t 32']
    )

    provir_fna=f'{checkv_dir}/proviruses.fna'
    vir_fna=f'{checkv_dir}/viruses.fna'
    combined_fna=f'{checkv_dir}/combined.fna'
    find_vir_cmd.extend(
        ['cat', provir_file, vir_file, combined_fna, '\n']
    )

    vs2_pass2_dir = f'{identify_dir}/vs2-pass2'
    if not os.path.exists(vs2_pass2_dir): os.makedirs(vs2_pass2_dir)
    find_vir_cmd.extend(
        ['virsorter', 'run', '--seqname-suffix-off', '--viral-gene-enrich-off',         '--provirus-off', '--prep-for-dramv', '-i', combined_fna, 
         '-w', vs2_pass2_dir, '--include-groups dsDNAphage,ssDNA',
         '--min-length 5000 --min-score 0.5 -j 28 all\n'
    )

    dramv_annotate_dir = f'{identify_dir}/dramv-annotate'
    if not os.path.exists(dramv_annotate_dir): os.makedirs(dramv_annotate_dir)
    vs2_pass2_fasta = f'{vs2_pass2_dir}/for-dramv/final-viral-combined-for-dramv.fa'
    vs2_pass2_tab = f'{vs2_pass2_dir}/for-dramv/viral-affi-contigs-for-dramv.tab'
    find_vir_cmd.extend(
        ['DRAM-v.py annotate', '-i', vs2_pass2_fasta, '-v', vs2_pass2_tab,  
         '-o', dramv_annotate_dir, '--skip_trnascan --threads 28 --min_contig_size 1000\n'
    )

    dramv_distill_dir = f'{identify_dir}/dramv-distill'
    dramv_annot = f'{dramv_annotate_dir}/annotations.tsv'
    find_vir_cmd.extend(
        ['DRAM-v.py distill', '-i', dramv_annot, '-o', dramv_distill_dir]  
    )

    curation_dir = f'{identify_dir}/curation'
    vir_score = f'{vs2-pass2_dir}/final-viral-score.tsv'
    curation_score = f'{curation_dir}/final-viral-score.tsv'
    contamination = f'{checkv_dir}/contamination.tsv'
    find_vir_cmd.extend(
        ["sed '1s/seqname/contig_id/'", vir_score, '>', curation_score, '\n']
    )
    


    



    


