#!/usr/bin/env python3
#coding=utf-8

import os
import sys
sys.path.append(sys.path[0]+'/crafts')
from crafts.dataqc import fastqc
from crafts.assembly import assembly
from crafts.identify import viridsop
from crafts.config import arguments
from crafts.votus import votus
from crafts.classify import classify
from crafts.quantify import align,quantify
from crafts.func_annot import geneAnnot
from crafts.compare import vCont

version='0.0.3'
args=arguments.setOpts(sys.argv[0],sys.argv[1],version)

if sys.argv[1]=='reads_qc':
    print('VirCraft data QC...\n')
    Reads=fastqc.QualCtrl(
        fq1=args.fq1,fq2=args.fq2,
        outdir=args.outdir,
        threads=args.threads
    )
    Reads.readqc(args.process)
    print('Reads quality control completed!!!')
elif sys.argv[1]=='assembly':
    print('VirCraft assembly...\n')
    Draft=assembly.Assembly(
        fq1=args.fq1,fq2=args.fq2,
        outdir=args.outdir,
        threads=args.threads
    )
    Draft.Assemble(
        process=args.process,
        cutoff=5000
    )
    print('Reads assembly completed!!!')
elif sys.argv[1]=='identify':
    print('Viral contig identification')
    VirSeq=viridsop.VirScan(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    VirSeq.Identify()
    print('Config identification completed!!!')
elif sys.argv[1]=='votus':
    print('Remove the redundancy')
    vOTUs=deRep.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup()
    print('vOTU cluster completed!!!')
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    Taxa=taxAnnot.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify()
    print('Contigs classification completed!!!')
elif sys.argv[1]=='quantify':
    print('Viral abundance and diversity analysis')
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
