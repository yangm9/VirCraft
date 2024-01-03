#!/usr/bin/env python3
#coding=utf-8
#author: yangming, yangm@idsse.ac.cn
#Created on Tus Jan 2 15:38:37 2024
import os
import sys
from crafts.config import arguments
from crafts.config import config

version='0.0.11'
parser=arguments.autOpts(sys.argv[0],version)
args=parser.parse_args()
if len(sys.argv)==1:
    parser.print_help()

samp_info = args.info
outdir = args.outdir

if not os.path.exists(outdir): os.makedirs(outdir)

CONFIG=config.VirCfg()
CONFIG.readSampInfo()
'''
