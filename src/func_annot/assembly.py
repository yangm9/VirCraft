#!/usr/bin/env python3

import sys,os
from ..config import setVari,conf
from ..process import cmdExec, general
from ..fastqc import mergeRead

def spades(fastqs: list, group: str, outdir: str):
    '''
    Assemble metagenome by SPAdes for single group.
    '''
    envs = setVari.selectENV('VirCraft')
    spades_cmd = [envs]
    spades_dir = f'{outdir}/01.assembly/{group}'
    general.mkdir(spades_dir)
    spades_cmd.extend(['spades.py', '--pe1-1', fastqs[0], '--pe1-2', fastqs[1], '--careful', '-t 30 -m 1300 -k 21,33,55,77,99,127', '-o', spades_dir])
    spades_sh = f'{outdir}/01.assembly/{group}_spades.sh'
    general.printSH(spades_sh, spades_cmd)
    results = cmdExec.execute(spades_cmd)
    return results

def filtFastA(grp: str, outdir: str, cutoff: int):
    '''
    Filter the fasta sequence by length (cutoff).
    '''
    wkdir = f'{outdir}/01.assembly/{grp}'
    scaffolds = f'{wkdir}/scaffolds.fasta'
    filt_fa_prifix = f'{wkdir}/scaffolds.filt'
    filt_cmd = ['SeqLenCutoff.pl', scaffolds, filt_fa_prifix, cutoff]
    filt_sh = f'{wkdir}/filt_scaffolds.sh'
    general.printSH(filt_sh, filt_cmd)
    results = cmdExec.execute(filt_cmd)
    return results

def statFastA(grp: str, outdir: str):
    wkdir = f'{outdir}/01.assembly/{grp}'
    contigs = f'{wkdir}/scaffolds.fasta'
    scaffolds = f'{wkdir}/scaffolds.fasta'
    stat_tab = f'{wkdir}/stat.tab'
    stat_cmd = ['assemb_stat.pl', contigs, scaffolds, f'>{stat_tab}\n']
    stat_sh = f'{wkdir}/stat_fasta.sh'
    general.printSH(stat_sh, stat_cmd)
    results = cmdExec.execute(stat_cmd)
    return results

def Assemble(config: str, outdir: str):
    mergeRead.MergeFastqByGroup(config,outdir)
    results=''
    for grp in mergeRead.groups:
        fastq_1 = f'{outdir}/01.assembly/0.fastq/{grp}_1.fq'
        fastq_2 = f'{outdir}/01.assembly/0.fastq/{grp}_2.fq'
        fastqs = [fastq_1, fastq_2]
        results += spades(fastqs, grp, outdir)
        results += filtFastA(grp, outdir, '2000')
        results += filtFastA(grp, outdir, '5000')
        results += filtFastA(grp, outdir, '10000')
        results += statFastA(grp, outdir)
    return results
