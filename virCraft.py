#!/usr/bin/env python
import sys
import OptionParser
from src.docs import checkOpts,getPara,readme,outPut
from src.assemble import assembly

usage=readme.description(sys.argv[0],'v0.0.1')
parser=OptionParser(usage)
parser.add_option(
    '-c','--configs',action='store',type='str',
    dest='config',metavar='STR',default=False,
    help='Cogfiguration file.'
)
parser.add_option(
    '-o','--output_dir',action='store',type='str',
    dest='output_dir',metavar='STR',default=False,
    help='Output direcortory.'
)
opts,args=parser.parse_args()

if sys.argv[1]=='assembly':
    print('VirCraft assembly')
    assembly.assemble(opts.config,opts.output_dir)
    exit(0)
    


