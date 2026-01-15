#!/usr/bin/env python3

import sys
import pandas as pd

# filt_checkv
# input: quality_summary.tsv, namely checkv result file
# methods: the vs2-vb-dvf-gn identified viral contigs with a length of ≥ 2 kb and .
# output: virus_quality_summary.filt.tsv and provir_quality_summary.filt.tsv
# return: viral ID list and proviral ID list
def filt_checkv(qual_summ_tsv: str):
    df = pd.read_csv(qual_summ_tsv, sep='\t')
    vctg_for_binning_filt_condition = 'contig_length >= 2000 and provirus != "Yes" and checkv_quality != "Complete" and checkv_quality != "Not-determined"'
    complete_vctg_filt_condition = f'provirus == "No" and checkv_quality == "Complete" and contig_length >= 2000'
    provirus_vctg_filt_condition = 'provirus == "Yes" and contig_length >= 2000 and checkv_quality != "Not-determined"'
    complete_vctg_list = df.query(complete_vctg_filt_condition)['Contig'].tolist()
    provirus_vctg_list = df.query(provirus_vctg_filt_condition)['Contig'].tolist()
    vctg_for_binning_list = df.query(vctg_for_binning_filt_condition)['Contig'].tolist()
    return complete_vctg_list, provirus_vctg_list, vctg_for_binning_list

# extracr sequences from FASTA file according to ID list
def extr_ctgs(id_list, in_fa, out_fa):
    IFA = open(in_fa)
    OFA = open(out_fa,'w')
    seq_id = ''
    curr_id = ''
    for line in IFA:
        line = line.strip()
        if line.startswith('>'):
            seq_id = line.lstrip('>')
            if seq_id in id_list:
                OFA.write(f'{line}\n')
                curr_id = seq_id
        else:
            if seq_id == curr_id:
                OFA.write(f'{line}\n')
    IFA.close()
    OFA.close()
    return 0

def filt_posi_ctgs(qual_summ_tsv: str, viral_filt_ctg_fna: str, viral_posi_ctg_fna: str):
    complete_vctg_list, provirus_vctg_list, vctg_for_binning_list = filt_checkv(qual_summ_tsv)
    virus_id_list = complete_vctg_list + provirus_vctg_list + vctg_for_binning_list
    extr_ctgs(virus_id_list, viral_filt_ctg_fna, viral_posi_ctg_fna)
    complete_vctg_fna = viral_posi_ctg_fna.replace('.fna', '_complete.fna')
    extr_ctgs(complete_vctg_list, viral_filt_ctg_fna, complete_vctg_fna)
    provirus_vctg_fna = viral_posi_ctg_fna.replace('.fna', '_provirus.fna')
    extr_ctgs(provirus_vctg_list, viral_filt_ctg_fna, provirus_vctg_fna)
    viral_binning_ctg_fna = viral_posi_ctg_fna.replace('.fna', '_for_binning.fna')
    extr_ctgs(vctg_for_binning_list, viral_filt_ctg_fna, viral_binning_ctg_fna)
    return 0

if __name__  == '__main__':
    if len(sys.argv) == 4:
        filt_posi_ctgs(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print(f'Usage: {sys.argv[0]} quality_summary.tsv viral_filt_ctg.fna viral_positive_ctg.fna')
