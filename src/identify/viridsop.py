
from ..config import setVari,conf
from ..process import cmdExec,general


def findVir(grp: str, outdir: str):
    envs = setVari.selectENV('VirCraft')
    find_vir_cmd = [envs]
    identify_dir = f'{outdir}/02.identify/{grp}'
    general.mkdir(identify_dir)
    vs2_pass1_dir = f'{identify_dir}/vs2-pass1'
    general.mkdir(vs2_pass1_dir)
    scaffolds = f'{outdir}/01.assembly/{grp}/scaffolds.filt.gt5000.fasta'
    find_vir_cmd.extend(
        ['virsorter', 'run', '--keep-original-seq', '-i', scaffolds, 
         '-d ', confDict['virsorter2DB'], '-v', vs2_pass1_dir, 
         '--include-groups dsDNAphage,ssDNA', 
         '--min-length 0 --min-score 0.5 -j 28 all\n']
    )

    vs2_pass1_fasta=f'{vs2_pass1_dir}/final-viral-combined.fa'
    checkv_dir=f'{identify_dir}/checkv'
    general.mkdir(checkv_dir)
    find_vir_cmd.extend(
        ['checkv', 'end_to_end', vs2_pass1_fasta, checkv_dir,
         '-d', confDict['CheckVDB'], '-t 32']
    )

    provir_fna=f'{checkv_dir}/proviruses.fna'
    vir_fna=f'{checkv_dir}/viruses.fna'
    combined_fna=f'{checkv_dir}/combined.fna'
    find_vir_cmd.extend(
        ['cat', provir_fna, vir_fna, combined_fna, '\n']
    )

    vs2_pass2_dir = f'{identify_dir}/vs2-pass2'
    general.mkdir(vs2_pass2_dir)
    find_vir_cmd.extend(
        ['virsorter', 'run', '--seqname-suffix-off', '--viral-gene-enrich-off',         '--provirus-off', '--prep-for-dramv', '-i', combined_fna, 
         '-w', vs2_pass2_dir, '--include-groups dsDNAphage,ssDNA',
         '--min-length 5000 --min-score 0.5 -j 28 all\n']
    )

    dramv_annotate_dir = f'{identify_dir}/dramv-annotate'
    general.mkdir(dramv_annotate_dir)
    vs2_pass2_fasta = f'{vs2_pass2_dir}/for-dramv/final-viral-combined-for-dramv.fa'
    vs2_pass2_tab = f'{vs2_pass2_dir}/for-dramv/viral-affi-contigs-for-dramv.tab'
    find_vir_cmd.extend(
        ['DRAM-v.py annotate', '-i', vs2_pass2_fasta, '-v', vs2_pass2_tab,  
         '-o', dramv_annotate_dir, '--skip_trnascan --threads 28 --min_contig_size 1000\n']
    )

    dramv_distill_dir = f'{identify_dir}/dramv-distill'
    general.mkdir(dramv_distill_dir)
    dramv_annot = f'{dramv_annotate_dir}/annotations.tsv'
    find_vir_cmd.extend(
        ['DRAM-v.py distill', '-i', dramv_annot, '-o', dramv_distill_dir, '\n']  
    )

    curation_dir = f'{identify_dir}/curation'
    general.mkdir(curation_dir)
    vir_score = f'{vs2_pass2_dir}/final-viral-score.tsv'
    curation_score = f'{curation_dir}/final-viral-score.tsv'
    contamination = f'{checkv_dir}/contamination.tsv'
    find_vir_cmd.extend(
        ["sed '1s/seqname/contig_id/'", vir_score, '>', curation_score, '\n']
    )
    
    cura_vs2_chkv = f'{curation_dir}/curation_vs2_checkv.tsv'
    find_vir_cmd.extend(
        ['linkTab.py', curation_score, contamination, 'left contig_id', cura_vs2_chkv, '\n']
    )
    
    find_vir_cmd.extend(
        ['vCurator.py', identify_dir, '\n']
    )
    combined_modi_fna = f'{checkv_dir}/combined_modi.fna'
    find_vir_cmd.extend(
        ["sed 's/_1 / /'", combined_fna, '>', combined_modi_fna, '\n']
    )

    contig_id_list = f'{curation_dir}/contigs_id.list'
    virus_posi_fna='{curation_dir}/virus_positive.fna'
    find_vir_cmd.extend(
        ['extrSeqByName.pl', contig_id_list, combined_modi_fna, virus_posi_fna, '\n']
    )
    find_vir_sh=f'{identify_dir}/find_vir.sh'
    general.printSH(find_vir_sh, find_vir_cmd)
    results = cmdExec.execute(find_vir_cmd)
    return results

def Identify(config: str, outdir: str):
    global confDict, sampDict
    groups, confDict, sampDict = conf.prepInfo(config)
    results = ''
    for grp in groups:
        results += f'{grp}: \n'
        results += findVir(grp, outdir)
    return results
