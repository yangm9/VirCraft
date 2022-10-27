#!/usr/bin/env python
import sys
from os import path
import argparse
from src.docs import readme
from src.assemble import assembly
from src.identify import viridsop
from src.process import general
from src.votus import deRep

usage = readme.description(sys.argv[0], VERSION)
parser = argparse.ArgumentParser(
    prog = 'VirCraft',
    usage = '%(prog)s <module> [opts] -o <output_file>',
    description = 'description: A flexible pipeline for metaviromic data analysis.'
)
parser.add_argument(
    '-c', '--configs', action='store', type=str,
    dest='config', metavar='STR', default=False,
    help='Cogfiguration file.'
)

parser.add_argument(
    '-t', '--data_type', action='store', type=str,
    dest='data_type', metavar='STR', default='FastQ',
    help='Input data type fastq orfasta.'
)

parser.add_argument(
    '-o', '--outdir', action='store', type=str,
    dest='outdir', metavar='STR', default=False,
    help='Output direcortory.'
)

args = parser.parse_args()

if args.outdir:
    outdir = path.abspath(args.outdir)
    general.mkdir(args.outdir)

if len(sys.argv) < 2:
    parser.print_help()
    exit(0)

if sys.argv[1] == 'assembly':
    print('VirCraft assembly')
    assembly.Assemble(args.config, outdir)
    exit(0)
elif sys.argv[1] == 'identify':
    print('Viral contig identification')
    viridsop.Identify(args.config, outdir)
    exit(0)
elif sys.argv[1] == 'votus':
    print('Remove the redundancy')
    deRep.RmDup(args.config, outdir)
    exit(0)
elif sys.argv[1] == 'classify':
    print('Viral contig classification')
    exit(0)
elif sys.argv[1] == 'diversity':
    print('Viral abundance and diversity')
    exit(0)
elif sys.argv[1] == 'function':
    print('Function annotation')
    exit(0)
elif sys.argv[1] == 'comparison':
    print('Viral ')
    exit(0)
elif sys.argv[1] == 'hosts':
    exit(0)
else:
    parser.error(f'{sys.argv[1]} is not a module of VirCraft!')
    exit(0)
