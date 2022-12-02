#!/usr/bin/env python3
#coding=utf-8

import os
import sys
sys.path.append(sys.path[0]+'/src')
from src.fastqc import readsQC
from src.assemble import assembly
from src.identify import viridsop
from src.config import argParser
from src.votus import deRep
from src.classify import taxAnnot
from src.quantify import align,coverage 
from src.func_annot import geneAnnot
from src.compare import vCont

version='0.0.3'
args=argParser.setOpts(sys.argv[0],version)

if sys.argv[1]=='reads_qc':
    print('VirCraft data QC...\n')
    Reads=readsQC.QualCtrl(fq1=args.fq1,fq2=args.fq2,outdir=args.outdir)
    Reads.readqc()
    print('Reads quality control completed!!!')
elif sys.argv[1]=='assembly':
    print('VirCraft assembly...\n')
    Draft=assembly.Assembly(args.fq1,args.fq2,args.outdir)
    Draft.Assemble()
    print('Reads assembly completed!!!')
elif sys.argv[1]=='identify':
    print('Viral contig identification')
    viridsop.Identify(opts.config,outdir)
elif sys.argv[1]=='votus':
    print('Remove the redundancy')
    deRep.RmDup(opts.config,outdir)
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    taxAnnot.demovir(opts.config,outdir)
elif sys.argv[1]=='quantify':
    print('Viral abundance and diversity')
    quantify.quantVir(opts.config,outdir)
elif sys.argv[1]=='func_annot':
    print('Function annotation')
    geneAnnot.funcAnno(opts.config,outdir)
elif sys.argv[1]=='comparison':
    print('Viral comparison')
    vCont.compSeq(opts.config,outdir)
elif sys.argv[1]=='host_prid':
    print('Host prediction')
else:
    ERROR=f'\nERROR: {sys.argv[1]} is not a module of VirCraft\n'
    print(ERROR)
    exit(0)
