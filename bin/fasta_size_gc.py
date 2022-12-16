#!/usr/bin/env python3
import sys
from Bio import SeqIO

def seqGcLen(seq:str):
    length=len(seq)
    gc_perc=round((seq.count('C')+seq.count('G'))/length*100)
    return length,gc_perc

def statFasta(fasta:str):
    print('SeqID\tGCcontent\tLength\n')
    for record in SeqIO.parse(fasta,'fasta'):
        length,gc_perc=seqGcLen(str(record.seq))
        length=str(length)
        gc_perc=str(gc_perc)
        print(f'{record.id}\t{length}\t{gc_perc}')
    return 0

if __name__=='__main__':statFasta(sys.argv[1])
