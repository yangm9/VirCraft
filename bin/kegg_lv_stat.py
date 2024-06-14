#!/usr/bin/env python3
import sys

def kegg_lv_stat(kegg_anno_f:str,lv_name:str):
    with open(kegg_anno_f) as f:
        gene_dict={}
        headers=f.readline().split('\t')
        col_idx=headers.index(lv_name)
        for line in f:
            items=line.strip().split('\t')
            gene_id=items[0]
            level=items[col_idx]
            gene_dict.setdefault(level,[])
            if not gene_id in gene_dict[level]:
                gene_dict[level].append(gene_id)
    for level,gene_id_list in gene_dict.items():
        set_length=len(set(gene_id_list))
        print(level+'\t'+str(set_length))
    return 0

if __name__=='__main__':
    kegg_lv_stat(sys.argv[1],sys.argv[2])
