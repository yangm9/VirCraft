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
from crafts.votus import uniqV
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
    ENV=installENV.ENV(
        outdir=args.outdir,
        threads=args.threads
    )
    ENV.Install(
        unrun=args.unrun,
        in_wall=args.in_wall
    )
    return 0

@logger.Log(level='INFO')
def setup_db(args):
    DB=deployDB.DB(
        outdir=args.outdir,
        threads=args.threads
    )
    DB.Deploy(
        unrun=args.unrun    
    )
    return 0

@logger.Log(level='INFO')
def reads_qc(args):
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
    return 0

@logger.Log(level='INFO')
def assemble(args):
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
    return 0

@logger.Log(level='INFO')
def identify(args):
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
    return 0

@logger.Log(level='INFO')
def votus(args):
    vOTUs=uniqV.VirRef(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    vOTUs.RmDup(
        args.cutoff,
        unrun=args.unrun,
        method=args.method
    )
    return 0

@logger.Log(level='INFO')
def classify(args):
    Taxa=vClassify.VirTaxa(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    Taxa.Classify(
        unrun=args.unrun
    )
    return 0

@logger.Log(level='INFO')
def compare(args):
    NWK=vCont.EnviComp(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads,
        orfprefix=args.orfprefix
    )
    NWK.CompSeq(
        unrun=args.unrun
    )
    return 0

@logger.Log(level='INFO')
def host_pred(args):
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
    AMG=callAMGs.AMGs(
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    AMG.annotAMGs(
        unrun=args.unrun
    )
    return 0

@logger.Log(level='INFO')
def vir_quant(args):
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
    return 0

@logger.Log(level='INFO')
def gene_quant(args):
    GeneQuant=geneQuantStat.GeneAbdStat(
        samp_info=args.samp_info,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads
    )
    GeneQuant.QuantStat(
        unrun=args.unrun
    )
    return 0

def main():
    version='0.0.14'
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
