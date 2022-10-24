#!/usr/bin/env python3
import sys,os,re
from ..config import setVari,conf
from ..process import cmdExec

def spades(fastq,outdir):
    conda_sh=
    envs=setVari.selectENV(conda_sh,VirCraft)
    spades_cmd=['spades.py','--pe1-1',fastq_1,'--pe1-2',fastq_2,'--careful','-t 30 -m 1300 -k 21,33,55,77,99,12','-o',group_name]
    print(' '.join(spades_cmd))
    results=cmdExec.execute(spades_cmd)
    return results

def Assemble(config: str,outdir: str):
    confDict,sampDict,groupDict=conf.prepInfo(config)
    
