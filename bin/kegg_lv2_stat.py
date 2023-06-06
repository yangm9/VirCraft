#!/usr/bin/env python3
import sys
with open(sys.argv[1]) as f:
    level2_sums={}
    f.readline()
    for line in f:
        columns=line.strip().split('\t')
        level2=columns[4]
        genes_annoted_interm=int(columns[1])
        if level2 in level2_sums:
            level2_sums[level2]+=genes_annoted_interm
        else:
            level2_sums[level2]=genes_annoted_interm
for level2,genes_annoted_interm_sum in level2_sums.items():
    print(level2+'\t'+str(genes_annoted_interm_sum))
