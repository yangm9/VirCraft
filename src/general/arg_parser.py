#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Fri Dec 18 18:17:58 2022
@author: yangming
'''

import argparse
#import src.version


def addFqArg(parser):
    parser.add_argument(
        '-1','--fastq1',action='store',type=str,
        dest='fq1',metavar='STR',default=False,
        required=True,help='FastQ file for read 1'
    )
    parser.add_argument(
        '-2','--fastq2',action='store',type=str,
        dest='fq2',metavar='STR',default=False,
        help='FastQ file for read 1'
    )
    return parser
def addSampFqArg(parser):
    parser.add_argument(
        '-q','--fastqs',action='store',type=str,
        dest='fqs',metavar='STR',default=False,
        help='Directory stored all FastQs'
    )
def addGlbArg(parser): #Add global arguments
    subparser.add_argument(
        '-t','--threads',action='store',type=str,
        dest='threads',metavar='INT',default=False,
        help='Num processes/threads to use\n'
    )
    subparser.add_argument(
        '-o','--outdir',action='store',type=str,
        dest='outdir',metavar='STR',default=False,
        required=True,help='Output direcortory\n'
    )
    return parser
#ver=src.version.__version__
parser=argparse.ArgumentParser(
    prog='VirCraft',
    description='VirCraft is an flexible pipeline for metaviromic data analysis.',
    usage=f'''
    virCraft.py <subcommands> [<options>] -o <outdir>
    subcommands: an optional functional module, including assembly, identify, votus, classify, quantify, func_annot and host_prid.
    options: options described below in the section of Options.
    outdir: output directory.
''',
    epilog = 'Text at the bottom of help'
)

subparsers=parser.add_subparsers(
    title='subcommands',
    description='valid subcommands',
    help=''
)
subparser=subparsers.add_parser('reads_qc',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('assembly',help='An assembly module could assemble the reads to contigs or scaffolds using MegaHit or SPAdes.')
subparser=addFqArg(subparser)
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('identify',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('votus',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('classify',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('quantify',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('func_annot',help='sub-command help')
subparser=addGlbArg(subparser)
subparser=subparsers.add_parser('host_prid',help='sub-command help')
subparser=addGlbArg(subparser)
args=parser.parse_args()
#args.func(args)
