#!/usr/bin/env python3
#coding=utf-8
#author: yangming, yangm@idsse.ac.cn
#Created on Fri Dec 18 18:17:58 2022

import os
import sys
sys.path.append(sys.path[0]+'/crafts')
from crafts.setup import installENV
from crafts.setup import deployDB
from crafts.data import fastqc
from crafts.assembly import assembly
from crafts.identify import viridsop
from crafts.identify import findV
from crafts.config import arguments
from crafts.votus import votus
from crafts.classify import classify
from crafts.classify import vCont
from crafts.quantify import virQuantStat
from crafts.quantify import geneQuantStat
from crafts.host import hosts
from crafts.func_annot import geneAnnot
from crafts.func_annot import callAMGs

version='0.0.11'
try:
    if len(sys.argv)<2:
        raise ValueError('''Insufficient parameters provided.
Please use -h or --help for assistance!''')
except ValueError as e:
    print(f'ERROR: {e}')
    exit(1)

args=arguments.setOpts(sys.argv[0],sys.argv[1],version)

if sys.argv[1]=='setup_env':
    print('Install the environments for VirCraft...\n')
    ENV=installENV.ENV(
        outdir=args.outdir,
        threads=args.threads
    )
    ENV.Install(
        unrun=args.unrun,
        in_wall=args.in_wall
    )
elif sys.argv[1]=='setup_db':
    print('VirCraft environments Done!!!\n')
    DB=deployDB.DB(
        outdir=args.outdir,
        threads=args.threads
    )
    DB.Deploy(
        unrun=args.unrun    
    )
    print('VirCraft environments Done!!!\n')
elif sys.argv[1]=='reads_qc':
    print('VirCraft data QC...\n')
    Reads=fastqc.QualCtrl(
        fq1=args.fq1,fq2=args.fq2,
        outdir=args.outdir,
        threads=args.threads
    )
    Reads.readqc(
        process=args.process,
        unrun=args.unrun,
        clear=args.clear
    )
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
        cutoff=args.cutoff,
        unrun=args.unrun
    )
    print('Reads assembly completed!!!')
elif sys.argv[1]=='identify':
    print('Viral contig identification')
    if args.sop=='viral-id-sop':
        VirSeq=viridsop.VirScan(
            fasta=args.fasta,
            outdir=args.outdir,
            threads=args.threads
        )
        VirSeq.Identify(
            unrun=args.unrun
        )
    else:
        VirSeq=findV.vIdentify(
            fasta=args.fasta,
            outdir=args.outdir,
            threads=args.threads
        )
        VirSeq.Identify(
            cutoff=args.cutoff,
            unrun=args.unrun
        )
    print('Config identification completed!!!')
elif sys.argv[1]=='votus':
    print('Remove the redundancy')
    vOTUs=votus.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup(
        args.cutoff,
        unrun=args.unrun
    )
    print('vOTU cluster completed!!!')
elif sys.argv[1]=='func_annot':
    #print('Function annotation')
    #VirGene=geneAnnot.GeneFunc(
    #    fasta=args.fasta,
    #    outdir=args.outdir,
    #    threads=args.threads
    #)
    #VirGene.FuncAnnot()
    print('Predicting AMGs')
    AMG=callAMGs.AMGs(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    AMG.annotAMGs(
        unrun=args.unrun
    )
    print('Function annotation completed!!!')
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    Taxa=classify.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify(
        unrun=args.unrun
    )
    print('Contigs classification completed!!!')
elif sys.argv[1]=='compare':
    print('Viral contig comparasion')
    NWK=vCont.EnviComp(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads,
        orfprefix=args.orfprefix
    )
    NWK.CompSeq(
        unrun=args.unrun
    )
    print('Contigs comparasion completed!!!')
elif sys.argv[1]=='vir_quant':
    print('Viral abundance and diversity analysis')
    VirQuant=virQuantStat.VirAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    VirQuant.QuantStat(
        args.taxa,
        args.checkv,
        unrun=args.unrun
    )
    print('Viral quantifications completed!!!')
elif sys.argv[1]=='gene_quant':
    print('gene abundance analysis')
    GeneQuant=geneQuantStat.GeneAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    GeneQuant.QuantStat(
        unrun=args.unrun
    )
    print('gene quantifications completed!!!')
elif sys.argv[1]=='host_pred':
    print('Host prediction')
    Hosts=hosts.VirHost(
        fasta=args.fasta,
        hostsdir=args.hostsdir,
        outdir=args.outdir,
        threads=args.threads
    )
    Hosts.PredHosts(
        gtdbtk=args.gtdbtkdir,
        taxa_anno=args.taxa,
        unrun=args.unrun
    )
    print('Viral-host relationship prediction done!!!')
else:
    ERROR=f'\nERROR: {sys.argv[1]} is not a module of VirCraft\n'
    print(ERROR)
    exit(0)
