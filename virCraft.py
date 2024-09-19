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
from crafts.taxa import vClassify
from crafts.taxa import vCont
from crafts.quantify import virQuantStat
from crafts.quantify import geneQuantStat
from crafts.host import hosts
from crafts.func_annot import geneAnnot
from crafts.func_annot import callAMGs
from crafts.general import logger


@logger.Log(level='INFO')
def setup_env(args):
    print('Installing all environments for VirCraft...')
    ENV=installENV.ENV(
        outdir=args.outdir,
        threads=args.threads
    )
    ENV.Install(
        unrun=args.unrun,
        in_wall=args.in_wall
    )
    print('VirCraft environments Done!!!')
    return 0

@logger.Log(level='INFO')
def setup_db(args):
    print('Deploying all databases for VirCraft...')
    DB=deployDB.DB(
        outdir=args.outdir,
        threads=args.threads
    )
    DB.Deploy(
        unrun=args.unrun    
    )
    print('VirCraft databases Done!!!')
    return 0

@logger.Log(level='INFO')
def reads_qc(args):
    print('VirCraft reads QC...')
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
    return 0

@logger.Log(level='INFO')
def assemble(args):
    print('VirCraft assembly...')
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
    return 0

@logger.Log(level='INFO')
def identify(args):
    print('Viral contig identification...')
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
    return 0

@logger.Log(level='INFO')
def votus(args):
    print('Remove the redundancy...')
    vOTUs=votus.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup(
        args.cutoff,
        unrun=args.unrun,
        method=args.method
    )
    print('\nvOTU cluster completed!!!')
    return 0

@logger.Log(level='INFO')
def classify(args):
    print('Viral contig classification...')
    Taxa=vClassify.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify(
        unrun=args.unrun
    )
    print('Contigs classification completed!!!')
    return 0

@logger.Log(level='INFO')
def compare(args):
    print('Viral contig comparasion...')
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
    return 0

@logger.Log(level='INFO')
def host_pred(args):
    print('Host prediction...')
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
    return 0

@logger.Log(level='INFO')
def func_annot(args):
    #print('Function annotation')
    #VirGene=geneAnnot.GeneFunc(
    #    fasta=args.fasta,
    #    outdir=args.outdir,
    #    threads=args.threads
    #)
    #VirGene.FuncAnnot()
    print('Predicting AMGs...')
    AMG=callAMGs.AMGs(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    AMG.annotAMGs(
        unrun=args.unrun
    )
    print('Function annotation completed!!!')
    return 0

@logger.Log(level='INFO')
def vir_quant(args):
    print('Viral abundance and diversity analysis...')
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
    return 0

@logger.Log(level='INFO')
def gene_quant(args):
    print('Gene abundance analysis...')
    GeneQuant=geneQuantStat.GeneAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    GeneQuant.QuantStat(
        unrun=args.unrun
    )
    print('Gene quantifications completed!!!')
    return 0

@logger.Log(level='INFO')
def main():
    version='0.0.13'
    parser=arguments.setOpts(sys.argv[0],version)
    args=parser.parse_args()
    if len(sys.argv)==1: 
        parser.print_help()
        exit(0)
    moduleDict={
        'setup_env':setup_env,
        'setup_db':setup_db,
        'reads_qc':reads_qc,
        'assemble':assemble,
        'identify':identify,
        'votus':votus,
        'classify':classify,
        'compare':compare,
        'vir_quant':vir_quant,
        'gene_quant':gene_quant,
        'func_annot':func_annot,
        'host_pred':host_pred
    }
    moduleDict[sys.argv[1]](args)
    return 0

if __name__=='__main__':
    main()
