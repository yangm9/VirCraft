#!/usr/bin/env python
#coding=utf-8
import sys
from os import path
from optparse import OptionParser #该模块已经不在开发维护，改为argparse
#sys.path.append(sys.path[0]+'/src')
from src.docs import readme
from src.assemble import assembly
from src.identify import viridsop
from src.general import general
from src.votus import deRep
from src.classify import taxAnnot
from src.quantify import quantVir 
from src.func_annot import geneAnnot
from src.compare import vCont

VERSION='v0.0.1'
usage=readme.description(sys.argv[0],VERSION)
parser=OptionParser(usage,version='%prog {VERSION}')
parser.add_option(
    '-c','--configs', action='store', type='str',
    dest='config',metavar='STR', default=False,
    help='Cogfiguration file.'
)
parser.add_option(
    '-t','--data_type', action='store', type='str',
    dest='data_type',metavar='STR', default='FastQ',
    help='Input data type fastq orfasta.'
)
parser.add_option(
    '-o','--outdir', action='store', type='str',
    dest='outdir',metavar='STR', default=False,
    help='Output direcortory.'
)
opts,args=parser.parse_args()

if opts.outdir:
    outdir=path.abspath(opts.outdir)
    general.mkdir(opts.outdir)

if len(sys.argv) < 2:
    parser.print_help()
    exit(0)

if sys.argv[1]=='assemble':
    print('VirCraft assembly')
    assembly.Assemble(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='identify':
    print('Viral contig identification')
    viridsop.Identify(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='votus':
    print('Remove the redundancy')
    deRep.RmDup(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    taxAnnot.demovir(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='quantify':
    print('Viral abundance and diversity')
    quantify.quantVir(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='func_annot':
    print('Function annotation')
    geneAnnot.funcAnno(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='comparison':
    print('Viral comparison')
    vCont.compSeq(opts.config,outdir)
    exit(0)
elif sys.argv[1]=='host_prid':
    print('Host prediction')
    exit(0)
else:
    parser.error(f'{sys.argv[1]} is not a module of VirCraft!')
    exit(0)
