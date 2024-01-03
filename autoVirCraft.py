#!/usr/bin/env python3
#coding=utf-8
#author: yangming, yangm@idsse.ac.cn
#Created on Tus Jan 2 15:38:37 2024
import os
import argparse
from crafts.config import config

parser = argparse.ArgumentParser(description='Process sample information')
parser.add_argument('--info', required=True, help='Path to sample_info.xls')
parser.add_argument('--outdir', required=True, help='Output directory for results')

args = parser.parse_args()

samp_info = args.info
outdir = args.outdir

if not os.path.exists(outdir): os.makedirs(outdir)

CONFIG=config.VirCfg()
CONFIG.readSampInfo()
