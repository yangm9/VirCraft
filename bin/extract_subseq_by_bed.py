#!/usr/bin/env pythn3
import sys
from Bio import SeqIO

def extrSeq(sequence, start, end):
    return sequence[start-1:end]

def getSubfaByBed(fasta_file,bed_file,subfa_file):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    #output_file=os.path.splitext(fasta_file)[0]
    SUBFA=open(output_file,"w")
    with open(bed_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            for fasta in fasta_sequences:
                if fasta.id == chrom:
                    sequence = str(fasta.seq)
                    subsequence = extrSeq(sequence, start, end)
                    SUBFA.write(">" + name + "\n")
                    SUBFA.write(subsequence + "\n")
                    break
    fasta_sequences.close()
    return 0

if __name__=='__main__':
    getSubfaByBed(sys.argv[1],sys.argv[2],sys.argv[3])
