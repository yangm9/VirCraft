#!/usr/bin/env python3
#coding=utf-8

import os
import sys
sys.path.append(sys.path[0]+'/src')
from src.assemble import assembly
from src.identify import viridsop
from src.config import arg_parser
from src.votus import deRep
from src.classify import taxAnnot
from src.quantify import align,coverage 
from src.func_annot import geneAnnot
from src.compare import vCont

version='0.0.3'
args=argParser.setOpts(sys.argv[0],version)

if sys.argv[1]=='assembly':
    print('VirCraft assembly')
    Draft=assembly.Assembly(opts.fq1,opts.fq2,opts.outdir)
    Draft.Assemble()
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
    ERROR=f'{sys.argv[1]} is not a module of VirCraft'
    parser.error(ERROR)
    exit(0)
