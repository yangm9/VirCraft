#!/usr/bin/env python3

import os
import sys
from ..config import setVari
from ..process import cmdExec, general
from ..fastqc.reads import Reads

class Assembly(Reads):
    '''
    '''
    envs = setVari.selectENV('VirCraft')
    def __init__(self, config, outdir):
        Reads.__init__(self, config, outdir)
        self.wkdir = f'{self.outdir}/01.assembly'
    def spades(self, fastqs: list, group: str):
        '''
        Assemble metagenome by SPAdes for single group.
        '''
        spades_cmd = [self.envs]
        wkdir = f'{self.wkdir}/{group}'
        general.mkdir(wkdir)
        spades_cmd.extend(
            ['spades.py', '--pe1-1', fastqs[0], '--pe1-2', fastqs[1],
            '--careful', '-t 30 -m 1300 -k 21,33,55,77,99,127', 
            '-o', wkdir]
        )
        spades_sh = f'{self.wkdir}/{group}_spades.sh'
        general.printSH(spades_sh, spades_cmd)
        results = cmdExec.execute(spades_cmd)
        return results
    def filtFastA(self, grp: str, cutoff: int):
        '''
        Filter the fasta sequence by length (cutoff).
        '''
        wkdir = f'{self.wkdir}/{grp}'
        scaffolds = f'{wkdir}/scaffolds.fasta'
        filt_fa_prifix = f'{wkdir}/scaffolds.filt'
        filt_cmd = ['SeqLenCutoff.pl', scaffolds, filt_fa_prifix, cutoff]
        filt_sh = f'{wkdir}/filt_scaffolds.sh'
        general.printSH(filt_sh, filt_cmd)
        results = cmdExec.execute(filt_cmd)
        return results
    def statFastA(self, grp: str):
        wkdir = f'{self.wkdir}/{grp}'
        contigs = f'{wkdir}/scaffolds.fasta'
        scaffolds = f'{wkdir}/scaffolds.fasta'
        stat_tab = f'{wkdir}/stat.tab'
        stat_cmd = ['assemb_stat.pl', contigs, scaffolds, f'>{stat_tab}\n']
        stat_sh = f'{wkdir}/stat_fasta.sh'
        general.printSH(stat_sh, stat_cmd)
        results = cmdExec.execute(stat_cmd)
        return results
    @property
    def Assemble(self):
        results=''
        for grp in self.groups:
            fastq_1 = f'{self.fq_dir}/{grp}_1.fq'
            fastq_2 = f'{self.fq_dir}/{grp}_2.fq'
            fastqs = [fastq_1, fastq_2]
            results += self.spades(fastqs, grp, outdir)
            results += filtFastA(grp, outdir, '2000')
            results += filtFastA(grp, outdir, '5000')
            results += filtFastA(grp, outdir, '10000')
            results += statFastA(grp, outdir)
        return results
