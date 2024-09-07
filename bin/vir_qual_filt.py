#!/usr/bin/env python3
import sys
import pandas as pd

def filt_checkv(checkv_d:str):
    '''
    input: full path of checkv directory
    methods: the vs2-vb-dvf identified viral contigs that are ≥ 5 kb or those that are circular and ≥ 1.5 kb were considered as high-confidence viral contigs.
    output: virus_quality_summary.filt.tsv and provir_quality_summary.filt.tsv
    return: viral ID list and proviral ID list
    '''
    virus_qual_filt_condition='provirus=="No" and (contig_length>=5000 or (contig_length>1500 and checkv_quality=="Complete"))'
    provir_qual_filt_condition='provirus=="Yes" and (proviral_length>=5000 or (proviral_length>1500 and checkv_quality=="Complete"))'
    qual_summ_tsv=f'{sys.argv[1]}/quality_summary.tsv'
    virus_qual_summ_tsv=qual_summ_tsv.replace('quality_summary.tsv','virus_quality_summary.filt.tsv')
    provir_qual_summ_tsv=qual_summ_tsv.replace('quality_summary.tsv','provir_quality_summary.filt.tsv')
    df=pd.read_csv(qual_summ_tsv,sep='\t')
    virus_df=df.query(virus_qual_filt_condition)
    provir_df=df.query(provir_qual_filt_condition)
    virus_id_list=virus_df['contig_id'].tolist()
    provir_id_list=provir_df['contig_id'].tolist()
    virus_df.to_csv(virus_qual_summ_tsv,index=False,sep='\t')
    provir_df.to_csv(provir_qual_summ_tsv,index=False,sep='\t')
    return virus_id_list,provir_id_list

def extr_ctgs(id_list,in_fa,out_fa):
    IFA=open(in_fa)
    OFA=open(out_fa,'w')
    seq_id=''
    curr_id=''
    for line in IFA:
        line=line.strip()
        if line.startswith('>'):
            seq_id=line.lstrip('>')
            if in_fa.endswith('proviruses.fna'):
                seq_id=seq_id.rsplit('_',1)[0]
            if seq_id in id_list:
                OFA.write(f'{line}\n')
                curr_id=seq_id
        else:
            if seq_id==curr_id:
                OFA.write(f'{line}\n')
    IFA.close()
    OFA.close()
    return 0

def filt_posi_ctgs(checkv_d:str):
    virus_id_list,provir_id_list=filt_checkv(checkv_d)
    virus_fna=f'{sys.argv[1]}/viruses.fna'
    virus_filt_fna=virus_fna.replace('.fna','.filt.fna')
    extr_ctgs(virus_id_list,virus_fna,virus_filt_fna)
    provir_fna=f'{sys.argv[1]}/proviruses.fna'
    provir_filt_fna=provir_fna.replace('.fna','.filt.fna')
    extr_ctgs(provir_id_list,provir_fna,provir_filt_fna)
    return 0

if __name__=='__main__':
   filt_posi_ctgs(sys.argv[1])
