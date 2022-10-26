#!/usr/bin/env python
import sys
from os import path
from optparse import OptionParser
from src.docs import readme
from src.assemble import assembly
from src.identify import viridsop
from src.process import general
from src.votus import deRep

usage = readme.description(sys.argv[0], 'v0.0.1')
parser = OptionParser(usage)
parser.add_option(
    '-c', '--configs', action='store', type='str',
    dest='config', metavar='STR', default=False,
    help='Cogfiguration file.'
)
parser.add_option(
    '-t', '--data_type', action='store', type='str',
    dest='data_type', metavar='STR', default='FastQ',
    help='Input data type fastq orfasta.'
)
parser.add_option(
    '-o', '--outdir', action='store', type='str',
    dest='outdir', metavar='STR', default=False,
    help='Output direcortory.'
)
opts, args = parser.parse_args()

general.mkdir(opts.outdir)
outdir = path.abspath(opts.outdir)

if sys.argv[1] == 'assembly':
    print('VirCraft assembly')
    assembly.Assemble(opts.config, outdir)
    exit(0)
elif sys.argv[1] == 'identify':
    print('Viral contig identification')
    viridsop.Identify(opts.config, outdir)
    exit(0)
elif sys.argv[1] == 'votu':
    print('Remove the redundancy')
    deRep(opts.config, outdir)
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
