#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Author:Yingli Zhou
# version 0.1.0 2021/1/18

import sys
from collections import defaultdict
from collections import Counter

"reference (https://doi.org/10.1016/j.cub.2017.10.040)"

VirusKeywordList = [
    "capsid","phage","terminase","base plate","baseplate", 
    "prohead","virion","virus","viral","tape measure",
    "tapemeasure neck","tail","head","bacteriophage","prophage",
    "portal","DNA packaging","T4","p22","holin"
]

def get_seq_dict(inputfile):
    "Read a FastA file，and return 'name':'sequence' dict."
    seq_dict={}
    with open(inputfile,'r') as f:
        for i in f.readlines():
            if i.startswith('>'):
                seqname=i.strip().strip('>')
                seqname=seqname.split(' ')[0]
                seq_dict[seqname]=''
            else:
                seq_dict[seqname]+=i.strip()
    return seq_dict

def get_contig_gene_num(protein_file):
    "Read a protein FastA file, and return 'name':'gene number' dict!"
    contig_gene_num_dict = {}
    with open(protein_file,'r') as f:
        for i in f.readlines():
            if i.startswith('>'):
                contig = i.strip().rsplit('_', 1)[0].strip('>')
                if contig not in contig_gene_num_dict.keys():
                    contig_gene_num_dict[contig] = 1
                else:
                    contig_gene_num_dict[contig] += 1
    return (contig_gene_num_dict)

def identify_phage(inputlist):
    "get virus specific gene list"
    for gene in inputlist:
        for i in VirusKeywordList:
            if gene.find(i)!=-1:
                virus_good_list.append(i)
    n=len(set(virus_good_list))
    if n>=2:
        return (True)

def get_positive_virus(eggfile,faa_f):
    positive_contig=[]
    virus_good_list=[]
    contig_dict=defaultdict(list)
    contig_protein_tax_dict=defaultdict(list)
    contig_protein_annotation_dict=defaultdict(list)
    virus_tax_dict = defaultdict(dict)
    with open(eggfile) as EGGF:
        for line in EGGF:
            if line.startswith('#query'): continue
            items=line.strip().split('\t')
            gene=items[0]
            family=items[4]
            dominant=items[17]
            annotation=items[-1]
            contig=gene.rsplit('_', 1)[0]
            if dominant=='Viruses':
                contig_protein_annotation_dict[contig].append('Viruses:'+annotation)
                contig_protein_tax_dict[contig].append(family)
            else:
                contig_protein_annotation_dict[contig].append(annotation)
    for contig in contig_dict.keys():
        if identify_phage(contig_dict[contig]):
            positive_contig.append(contig)
    contig_gene_num=get_contig_gene_num(faa_f)
    for contig in contig_protein_annotation_dict.keys():
        known_num=0
        for annotation in contig_protein_annotation_dict[contig]:
            if not annotation.startswith('Viruses:'):
                if annotation:
                    if annotation.find("unknown") != -1 or annotation.find("hypothetical protein") != -1:
                        continue
                    else:
                        known_num+=1
            if annotation.find('ribosomal protein') != -1:
                known_num = contig_gene_num[contig]
                break
        print(float(known_num)/contig_gene_num[contig])
        if float(known_num)/contig_gene_num[contig] > 0.3:
            if contig in positive_contig:
                positive_contig.remove(contig)
        else:
            if contig not in positive_contig:
                positive_contig.append(contig)
    return positive_contig

def main(
    fa_f,
    faa_f,
    eggout_f='eggout.emapper.annotations',
    out_f='virus_keywords_positive.fa'
):
    result=[]
    positive_virus=get_positive_virus(eggout_f,faa_f)
    seq_dict=get_seq_dict(fa_f)
    for contig in positive_virus:
        result.append(">%s\n%s\n" % (contig,seq_dict[contig]))
    with open(out_f,'w') as f:
        f.writelines(result)

if __name__ == '__main__':
    if len(sys.argv)>=3:
        main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        print('python '+sys.argv[0]+' <virus_contigs.fa> <virus_contigs.faa> <eggout.emapper.annotations> <output.fa>')
