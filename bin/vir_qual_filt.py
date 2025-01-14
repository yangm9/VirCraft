#!/usr/bin/env python3
import sys
import pandas as pd

def filt_checkv(qual_summ_tsv: str):
    '''
    input: quality_summary.tsv, namely checkv result file
    methods: the vs2-vb-dvf-gm identified viral contigs that are ≥ 5 kb or those that are circular and ≥ 1.5 kb were considered as high-confidence viral contigs.
    output: virus_quality_summary.filt.tsv and provir_quality_summary.filt.tsv
    return: viral ID list and proviral ID list
    '''
    qual_filt_condition = 'contig_length >= 5000 or (contig_length >= 1500 and checkv_quality == "Complete")'
    filt_qual_summ_tsv = qual_summ_tsv.replace('quality_summary.tsv', 'quality_summary.filt.tsv')
    df = pd.read_csv(qual_summ_tsv, sep='\t')
    filt_df = df.query(qual_filt_condition)
    filt_id_list = filt_df['contig_id'].tolist()
    filt_df.to_csv(filt_qual_summ_tsv, index = False, sep = '\t')
    return filt_id_list

def extr_ctgs(id_list,in_fa,out_fa):
    IFA = open(in_fa)
    OFA = open(out_fa,'w')
    seq_id = ''
    curr_id = ''
    for line in IFA:
        line=line.strip()
        if line.startswith('>'):
            seq_id=line.lstrip('>')
            if seq_id in id_list:
                OFA.write(f'{line}\n')
                curr_id=seq_id
        else:
            if seq_id==curr_id:
                OFA.write(f'{line}\n')
    IFA.close()
    OFA.close()
    return 0

def filt_posi_ctgs(qual_summ_tsv: str, viral_filt_ctg_fna: str, viral_posi_ctg_fna: str):
    virus_id_list = filt_checkv(qual_summ_tsv)
    extr_ctgs(virus_id_list, viral_filt_ctg_fna, viral_posi_ctg_fna)
    return 0

if __name__  == '__main__':
    if len(sys.argv) == 4:
        filt_posi_ctgs(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print(f'Usage: {sys.argv[0]} quality_summary.tsv viral_filt_ctg.fna viral_positive_ctg.fna')
