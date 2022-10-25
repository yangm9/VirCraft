#!/usr/bin/env python
import sys
from optparse import OptionParser
from src.docs import readme
from src.assemble import assembly
from src.identify import viridsop

usage=readme.description(sys.argv[0],'v0.0.1')
parser=OptionParser(usage)
parser.add_option(
    '-c','--configs',action='store',type='str',
    dest='config',metavar='STR',default=False,
    help='Cogfiguration file.'
)
parser.add_option(
    '-t','--data_type',action='store',type='str',
    dest='data_type',metavar='STR',default='FastQ',
    help='Input data type fastq orfasta.'
)
parser.add_option(
    '-o','--output_dir',action='store',type='str',
    dest='output_dir',metavar='STR',default=False,
    help='Output direcortory.'
)
opts,args=parser.parse_args()

if sys.argv[1]=='assembly':
    print('VirCraft assembly')
    assembly.Assemble(opts.config,opts.output_dir)
    exit(0)
elif sys.argv[1]=='identify':
    print('Viral contig identification')
    viridsop.findVir(opts.config,opts.output_dir)
    exit(0)
elif sys.argv[1]=='classify':
    print('Viral contig classification')
    exit(0)
elif sys.argv[1]=='diversity':
    print('Viral abundance and diversity')
    exit(0)
elif sys.argv[1]=='function':
    print('Function annotation')
    exit(0)
elif sys.argv[1]=='comparison':
    print('Viral ')
    exit(0)
elif sys.argv[1]=='hosts':
    exit(0)
