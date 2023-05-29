#!/usr/bin/env python3
import os
import sys
import pandas as pd
import linkTab

def get_fasta_seq_names(votu_f):
    seq_names = []
    with open(votu_f, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_names.append(line[1:])
    df = pd.DataFrame(seq_names, columns=['Sequence_ID'])
    return df

def complete_contigs(votu_f,blast_taxa_f):
    ctg_df=get_fasta_seq_names(votu_f)
    name=os.path.splitext(os.path.basename(votu_f))[0]
    outdir=os.path.dirname(blast_taxa_f)
    votu_list=f'{outdir}/temp_votu_id.list'
    ctg_df.to_csv(votu_list,sep='\t',index=False)
    out_f=f'{outdir}/{name}.taxa.txt'
    linkTab.merge(votu_list,blast_taxa_f,'left','Sequence_ID',out_f)
    return 0

if __name__=='__main__':
    complete_contigs(sys.argv[1],sys.argv[2])


