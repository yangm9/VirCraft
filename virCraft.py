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
#try:
#    if len(sys.argv)<2:
#        raise ValueError('''Insufficient parameters provided.
#Please use "virCraft.py -h" or "virCraft.py --help" for assistance!''')
#except ValueError as e:
#    print(f'ERROR: {e}')
#    exit(1)
parser=arguments.setOpts(sys.argv[0],version)
args=parser.parse_args()
if len(sys.argv)==1: 
    parser.print_help()
    exit(0)

if sys.argv[1]=='setup_env':
    print('Installing all environments for VirCraft...\n')
    ENV=installENV.ENV(
        outdir=args.outdir,
        threads=args.threads
    )
    ENV.Install(
        unrun=args.unrun,
        in_wall=args.in_wall
    )
    print('\nVirCraft environments Done!!!')
elif sys.argv[1]=='setup_db':
    print('Deploying all databases for VirCraft...\n')
    DB=deployDB.DB(
        outdir=args.outdir,
        threads=args.threads
    )
    DB.Deploy(
        unrun=args.unrun    
    )
    print('\nVirCraft databases Done!!!')
elif sys.argv[1]=='reads_qc':
    print('VirCraft reads QC...\n')
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
    print('\nReads quality control completed!!!')
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
    print('\nReads assembly completed!!!')
elif sys.argv[1]=='identify':
    print('Viral contig identification...\n')
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
    print('\nConfig identification completed!!!')
elif sys.argv[1]=='votus':
    print('Remove the redundancy...\n')
    vOTUs=votus.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup(
        args.cutoff,
        unrun=args.unrun
    )
    print('\nvOTU cluster completed!!!')
elif sys.argv[1]=='func_annot':
    #print('Function annotation')
    #VirGene=geneAnnot.GeneFunc(
    #    fasta=args.fasta,
    #    outdir=args.outdir,
    #    threads=args.threads
    #)
    #VirGene.FuncAnnot()
    print('Predicting AMGs...\n')
    AMG=callAMGs.AMGs(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    AMG.annotAMGs(
        unrun=args.unrun
    )
    print('\nFunction annotation completed!!!')
elif sys.argv[1]=='classify':
    print('Viral contig classification...\n')
    Taxa=classify.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify(
        unrun=args.unrun
    )
    print('\nContigs classification completed!!!')
elif sys.argv[1]=='compare':
    print('Viral contig comparasion...\n')
    NWK=vCont.EnviComp(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads,
        orfprefix=args.orfprefix
    )
    NWK.CompSeq(
        unrun=args.unrun
    )
    print('\nContigs comparasion completed!!!')
elif sys.argv[1]=='vir_quant':
    print('Viral abundance and diversity analysis...\n')
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
    print('\nViral quantifications completed!!!')
elif sys.argv[1]=='gene_quant':
    print('Gene abundance analysis...\n')
    GeneQuant=geneQuantStat.GeneAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    GeneQuant.QuantStat(
        unrun=args.unrun
    )
    print('\nGene quantifications completed!!!')
elif sys.argv[1]=='host_pred':
    print('Host prediction...\n')
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
    print('\nViral-host relationship prediction done!!!')
else: 
    exit(0)
