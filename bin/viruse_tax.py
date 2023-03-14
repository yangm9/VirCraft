#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author:Yingli Zhou
# modified by: yangming@idsse.ac.cn
# version 0.1.0 2022/1/4
# version 0.1.1 2023/3/13

import sys
from collections import defaultdict
from collections import Counter

contig_dict = defaultdict(list)
contig_family = defaultdict(list)
contig_class = defaultdict(list)

## get the dictionary of contig:[gene1,gene2...]
# with open('/data_backup/zhouyl/Project/1.Mariana/3.virus/2.Annotation/5.checkV/tmp/CD_viral.faa','r') as f:
with open(sys.argv[1],'r') as f:
    for i in f.readlines():
        if i.startswith('>'):
            gene = i.strip().strip('>')
            contig = gene.rsplit('_',1)[0]
            contig_dict[contig].append(gene)

## get the dictionary of contig:[class,family]
# with open('CD_gene_taxonomy.txt','r') as f:
with open(sys.argv[2],'r') as f:
    lastgene=''
    for i in f.readlines()[1:]:
        ipart = i.strip().split('\t')
        thisgene=ipart[0]
        contig = ipart[0].rsplit('_',1)[0]
        taxonomy = ipart[-1]
        class_level = taxonomy.split(';')[3]
        family = taxonomy.split(';')[4]
        if thisgene != lastgene:
            contig_class[contig].append(class_level)
            contig_family[contig].append(family)
        lastgene=thisgene

for contig in contig_dict.keys():
    gene_num = len(contig_dict[contig])
    class_rate = 0
    if contig in contig_family.keys():
        class_static = Counter(contig_class[contig])
        for class_level in class_static.keys():
            class_rate = class_static[class_level]/float(gene_num)
            if class_rate>=0.5:
                line= contig + '\t' + class_level
                break
        if class_rate <0.5:
            line = contig + '\tUnassigned'
    else:
        line = contig + '\tUnassigned'

    if contig in contig_family.keys():
        family_static = Counter(contig_family[contig])
        for family in family_static.keys():
            family_rate = family_static[family]/float(gene_num)
            if family_rate>=0.5:
                line = line + '\t' + family
                break
        if family_rate <0.5:
            line = line +'\tUnassigned'
    else:
        line = line + '\tUnassigned'

    print(line)

