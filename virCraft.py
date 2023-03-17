#!/usr/bin/env python3
#coding=utf-8

import os
import sys
sys.path.append(sys.path[0]+'/crafts')
from crafts.data import fastqc
from crafts.assembly import assembly
from crafts.identify import viridsop
from crafts.config import arguments
from crafts.votus import votus
from crafts.classify import classify
from crafts.classify import vCont
from crafts.quantify import virQuantStat
from crafts.host import hosts
from crafts.func_annot import geneAnnot

version='0.0.6'
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
        cutoff=args.cutoff
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
    vOTUs=votus.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup()
    print('vOTU cluster completed!!!')
elif sys.argv[1]=='vir_quant':
    print('Viral abundance and diversity analysis')
    VirQuant=virQuantStat.VirAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    VirQuant.QuantStat(args.taxa,args.checkv)
    print('Viral quantifications completed!!!')
elif sys.argv[1]=='func_annot':
    print('Function annotation')
    VirGene=geneAnnot.GeneFunc(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    VirGene.FuncAnnot()
    print('Function annotation completed!!!')
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    Taxa=classify.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify()
    print('Contigs classification completed!!!')
elif sys.argv[1]=='compare':
    print('Viral contig comparasion')
    NWK=vCont.EnviComp(
        orfs=args.fasta,
        outdir=args.outdir,
        threads=args.threads,
        orfprefix=args.orfprefix
    ) #network
    NWK.vContact()
    print('Contigs comparasion completed!!!')
elif sys.argv[1]=='host_prid':
    print('Host prediction')
    Hosts=hosts.VirHost(
        fasta=args.fasta,
        hostsdir=args.hostsdir,
        outdir=args.outdir,
        threads=args.threads
    )
    Hosts.PredHosts()
    print('Viral-host relationship prediction done!!!')
else:
    ERROR=f'\nERROR: {sys.argv[1]} is not a module of VirCraft\n'
    print(ERROR)
    exit(0)
