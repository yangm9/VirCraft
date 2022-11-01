#!/usr/bin/env python3
import sys

with open(sys.argv[1],'r') as f:
    print('protein_id,contig_id,keywords')
    for i in f.readlines():
        if i.startswith('>'):
            protein_id=i.strip().strip('>')
            contig_index=i.rfind('_')
            contig_id=i[1:contig_index]
            print("%s,%s,None_provided" % (protein_id,contig_id))
