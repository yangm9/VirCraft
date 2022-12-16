#!/usr/bin/env python3
from Bio import SeqIO

def gc_content(seq):
    return round((seq.count('C')+seq.count('G'))/len(seq)*100

def processFasta(fas):
    print("SeqID\tGCcontent\tLength\n")
    for record in SeqIO.parse(fas, "fasta"):
        print(f'{record.id\tgetGC(str(record.seq))\tstr(len(record.seq))))
