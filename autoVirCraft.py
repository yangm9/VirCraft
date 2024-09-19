#!/usr/bin/env python3
#coding=utf-8
#author: yangming, yangm@idsse.ac.cn
#Created on Tus Jan 2 15:38:37 2024
import os
import sys
from crafts.config import arguments
from crafts.config import config

version='0.0.14'
parser=arguments.autOpts(sys.argv[0],version)
args=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()

if not os.path.exists(outdir): os.makedirs(outdir)
if sys.argv[1]=='exec':
    samp_info=args.samp_info
    fasta=args.fasta
    threads=args.threads
    outdir=args.outdir
    CONFIG=config.VirCfg()
    groups,sampDict=CONFIG.readSampInfo(samp_info)
    threads_per_samp=threads/len(groups)
    reads_qc_dir=f'{outdir}/01.reads_qc'
    for samp in sampDict.keys():
        if '1' in steps:
            fq1,fq2=sampDict[samp][1].split(',')
            samp_outdir=f'{reads_qc_dir}/{samp}'
            cmd=['virCraft.py','reads_qc','-1',fq1,'-2',fq2,
                '-t',threads_per_samp,'-p fuc','-o',samp_outdir,'\n']
            shell=f'{reads_qc_dir}/{samp}_reads_qc.sh'
            utils.printSH(shell,cmd)
        else:

        if '2' in steps:

