#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Fri Dec 18 18:17:58 2022
@author: yangming
'''

import argparse
#import src.version

def setOpts(name:str,subcmds:str,version:str):
    #ver=src.version.__version__
    parser=argparse.ArgumentParser(
        prog=name,
        description='VirCraft is an flexible pipeline for metaviromic data analysis.',
        usage=f'''
        {name} {subcmds} [<options>] -o <outdir>
        subcommands: an optional functional module, including assembly, identify, votus, classify, quantify, func_annot and host_prid.
        options: options described below in the section of Options.
        outdir: output directory.
    ''',
        epilog='Text at the bottom of help'
    )

    #Create subcommands objects
    subparsers=parser.add_subparsers(
        title='subcommands',
        description='valid subcommands',
        help=''
    )
    
    subpsr=subparsers.add_parser(
        'reads_qc',
        help='Pair-end FastQ reads qualitiy control.'
    )
    subpsr=addPairFqArg(subpsr)
    subpsr=addGlbArg(subpsr)
    subpsr=addProcArg(subpsr,'fuc')
    
    subpsr=subparsers.add_parser(
        'assembly',
        help='Assemble the reads to contigs or scaffolds using MegaHit or SPAdes'
    )
    subpsr=addPairFqArg(subpsr)
    subpsr=addGlbArg(subpsr)
    subpsr=addProcArg(subpsr,'sm')

    subpsr=subparsers.add_parser(
        'identify',
        help='identify the viral contigs from a assembly fasta, using vir-id-sop'
    )
    subpsr=addFaArg(subpsr)
    subpsr=addGlbArg(subpsr)
    
    subpsr=subparsers.add_parser(
        'votus',
        help='construct the non-redundant virus operational taxonomic unit (vOTU) reference'
    )
    subpsr=addFaArg(subpsr)
    subpsr=addGlbArg(subpsr)
    
    subpsr=subparsers.add_parser(
        'classify',
        help='sub-command help'
    )
    subpsr=addFaArg(subpsr)
    subpsr=addGlbArg(subpsr)
    
    subpsr=subparsers.add_parser(
        'quantify',
        help='sub-command help'
    )
    subpsr=addAllFqArg(subpsr)
    subpsr=addFaArg(subpsr)
    subpsr=addGlbArg(subpsr)
    subpsr.add_argument(
        '-s','--sampinfo',action='store',type=str,
        dest='samp_info',metavar='STR',default=False,
        required=True,help='Taxonomic annotation file'
    )
    subpsr.add_argument(
        '-x','--taxa',action='store',type=str,
        dest='taxa',metavar='STR',default=False,
        help='Taxonomic annotation file'
    )
    
    subpsr=subparsers.add_parser(
        'func_annot',
        help='sub-command help'
    )
    subpsr=addFaArg(subpsr)
    subpsr=addGlbArg(subpsr)
    
    subpsr=subparsers.add_parser(
        'host_prid',
        help='sub-command help'
    )
    subpsr=addGlbArg(subpsr)
    
    args=parser.parse_args()
    return args

def addGlbArg(psr): #Add global arguments
    psr.add_argument(
        '-t','--threads',action='store',type=str,
        dest='threads',metavar='INT',default=8,
        help='Num processes/threads to use\n'
    )
    psr.add_argument(
        '-o','--outdir',action='store',type=str,
        dest='outdir',metavar='STR',default=False,
        required=True,help='Output direcortory\n'
    )
    return psr

def addProcArg(psr,dflt:str):
    psr.add_argument(
        '-p','--process',action='store',type=str,
        dest='process',metavar='STR',default=dflt,
        help='Select the analysis process for a certain module\n'
    )
    return psr

def addPairFqArg(psr):
    psr.add_argument(
        '-1','--fastq1',action='store',type=str,
        dest='fq1',metavar='STR',default=False,
        required=True,help='FastQ file for read 1'
    )
    psr.add_argument(
        '-2','--fastq2',action='store',type=str,
        dest='fq2',metavar='STR',default=False,
        help='FastQ file for read 2'
    )
    return psr

def addAllFqArg(psr):
    psr.add_argument(
        '-q','--fastqs',action='store',type=str,
        dest='fastqs',metavar='STR',default=False,
        help='All clean FastQs, i.e. path/*.fastq'
    )
    return psr

def addFaArg(psr):
    psr.add_argument(
        '-a','--fasta',action='store',type=str,
        dest='fasta',metavar='STR',default=False,
    help='a Fasta file'
    )
    return psr