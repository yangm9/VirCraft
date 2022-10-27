#!/usr/bin/env python3

import sys,os
import shutil
from ..config import setVari, conf
from ..process import cmdExec, general

def demovir(config: str, outdir: str):
    '''
    Classify the virus contig by Demovir software for a certain single group.
    '''
    groups, confDict, sampDict = conf.prepInfo(config)
    envs = setVari.selectENV('VirCraft')
    demovir_cmd = [envs]
    wkdir = f'{outdir}/04.classify'
    general.mkdir(wkdir)
    TrEMBL_viral_taxa = f'{confDict["DemovirDB"]}/TrEMBL_viral_taxa.RDS'
    TrEMBL_viral_taxa_lnk = f'{wkdir}/TrEMBL_viral_taxa.RDS'
    demovir_cmd.extend(['ln -s', TrEMBL_viral_taxa, TrEMBL_viral_taxa_lnk, '\n'])
    uniprot_trembl_viral = f'{confDict["DemovirDB"]}/uniprot_trembl.viral.udb'
    uniprot_trembl_viral_lnk = f'{wkdir}/uniprot_trembl.viral.udb'
    demovir_cmd.extend(['ln -s', uniprot_trembl_viral, uniprot_trembl_viral_lnk, '\n'])
    demovir = f'{sys.path[0]}/bin/demovir.*'
    demovir_cmd.extend(['cp', demovir, wkdir, '\n'])
    votus = f'{outdir}/03.vOTUs/merged_virus_positive_nodup.fa'
    demovir = f'{wkdir}/demovir.sh'
    demovir_cmd.extend([demovir, votus, '32'])
    demovir_sh = f'{wkdir}/classify_by_demovir.sh'
    general.printSH(demovir_sh, demovir_cmd)
    results = cmdExec.execute(demovir_cmd)
    return results
