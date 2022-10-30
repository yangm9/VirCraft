#!/usr/bin/env python3

import sys,os
from ..config import setVari,conf
from ..process import cmdExec, general
from ..fastqc import mergeRead

def mkIdx(fasta: str):
    envs = setVari.selectENV('VirCraft')
    idx_cmd = [envs]
    fasta_lnk = f'{wkdir}/merged_virus_positive_nodup.fa'
    idx_cmd.extend(['ln -s',fasta, fasta_lnki, '\n'])
    bwa_idx_sh = f'{wkdir}/bwa_index.sh'
    idx_cmd.extend(['bwa index -a bwtsw', fasta_lnk, '-p', bwa_idx, '\n'])
    general.printSH(bwa_idx_sh, idx_cmd)
    results = cmdExec.execute(spades_cmd)
    return results

def alnReads(samp: str):
    '''
    Align the reads to the vOTUs for each sample.
    '''
    envs = setVari.selectENV('VirCraft')
    bwa_cmd = [envs]
    raw_bam = f'{wkdir}/{samp}.raw.bam'
    sort_bam = f'{wkdir}/{samp}.sort.bam'
    fastqs = sampDict[samp][2].split(',')
    bwa_cmd.extend(
        ['bwa mem', '-t 28', bwa_idx, fastqs[0], fastqs[1], 
         '|samtools view', '-o', raw_bam, '-@ 28 -b -S\n', 
         'samtools sort', raw_bam, '-o', sort_bam, '-@ 28\n']
    )
    bwa_sh = f'{outdir}/0.bwa/{samp}_bwa.sh'
    general.printSH(bam_sh, bam_cmd)
    results = cmdExec.execute(bam_cmd)
    return results

def AlignBySamp(config: str, outdir: str):
    global sampDict, bwa_idx, wkdir
    groups, confDict, sampDict = conf.prepInfo(config)
    bwa_idx = f'{wkdir}/merged_virus_positive_nodup'
    wkdir = f'{outdir}/05.abundance/1.bwa'
    general.mkdir(wkdir)
    fasta = f'{outdir}/03.vOTUs/merged_virus_positive_nodup.fa'
    results = mkIdx(fasta)
    for samp in sampDict.keys():
        results += alnReads(samp) 
    return results
