#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on Fri Dec 18 18:17:58 2022
@author: yangming
'''

import argparse

def autOpts(name: str, version: str):
    parser = argparse.ArgumentParser(
        prog = name,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        usage = name,
        description = f'''autoVirCraft {version} -- A automatic program for VirCraft
Usage: {name} <subcommand> [options]
'''
    )
    subparsers=parser.add_subparsers(
        help=''
    )
    subpsr = subparsers.add_parser(
        'exec',
        help='Execute the VirCraft pipeline automatically'
    )
    subpsr = addSampArg(subpsr)
    subpsr = addFaArg(subpsr,'ctg')
    subpsr.add_argument(
        '-p', '--steps', action='store', type=str, dest='steps', metavar='STR', default='123456789', required=False,
        help='the steps selected for viromic analysis, [default=123456789]'
    )
    
    subpsr = subparsers.add_parser(
        'setup',
        help='Setup the conda environments and databases for VirCraft automatically'
    )
    
    subpsr = subparsers.add_parser(
        'clear',
        help='Remove all useless files from the results to save storage space'
    )
    return parser

def setOpts(name: str, version: str):
    parser = argparse.ArgumentParser(
        prog = name,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        #usage=argparse.SUPPRESS,
        usage = name,
        description = f'''VirCraft {version} -- A flexible pipeline for metaviromic data analysis

Usage: {name} <subcommand> [options] -o <outdir>
    subcommand: an optional functional module, including reads_qc, assemble, identify, votus, binning, classify, compare, host_pred, func_annot, vir_quant and gene_quant.
    options: options described below in the section of Options.
    outdir: output folder.
''',
        epilog = 'Text at the bottom of help'
    )

    #Create subcommands objects
    subparsers = parser.add_subparsers(
        help=''
    )
#----------------------setup_env-----------------------
    subpsr = subparsers.add_parser(
        'setup_env',
        help='Install the conda environments required by VirCraft'
    )
    subpsr = addGlbArg(subpsr)
    subpsr.add_argument(
        '-w', '--in-wall', action='store_true', dest='in_wall', default=False, required=False,
        help='Is the network of your server within the Great Wall? [default=False]'
    )

#----------------------setup_db-----------------------
    subpsr = subparsers.add_parser(
        'setup_db',
        help='Deploy the bioinformatic databases requireed by VirCraft'
    )
    subpsr = addGlbArg(subpsr)

#----------------------reads_qc-----------------------
    subpsr = subparsers.add_parser(
        'reads_qc',
        help='Pair-end FastQ reads qualitiy control via fastp, fastuniq and/or decontamination'
    )
    subpsr = addPairFqArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr = addProcArg(subpsr, 'fuc')
    
#----------------------assemble-----------------------
    subpsr = subparsers.add_parser(
        'assemble',
        help='Assemble the pair-end reads to contigs or scaffolds using MegaHit and/or SPAdes'
    )
    subpsr = addPairFqArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr = addProcArg(subpsr, 'ms')
    subpsr = addCutoffArg(subpsr)

#----------------------identify-----------------------
    subpsr = subparsers.add_parser(
        'identify',
        help='Identify the viral contigs from an assembly fasta, using vs2-vb-dvf-gm or vir-id-sop workflow'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr = addCutoffArg(subpsr)
    subpsr.add_argument(
        '-w', '--sop', action='store', type=str, dest='sop', metavar='STR', default='vs2-vb-dvf-gm', required=False,
        help='The sop/pipeline for viral contigs identification, including "viral-id-sop" and "vs2-vb-dvf-gm". [default=vs2-vb-dvf-gm]'
    )
    subpsr.add_argument(
        '-f', '--filter_mode', action='store', type=str, dest='mode', metavar='STR', default='permissive', required=False,
        help='Filter mode for viral contigs, including "permissive" and "strict". [default=permissive]'
    )
 
#-----------------------votus-------------------------
    subpsr = subparsers.add_parser(
        'votus',
        help='Construct the non-redundant virus operational taxonomic unit (vOTU) reference'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr = addCutoffArg(subpsr)
    subpsr.add_argument(
        '-m', '--method', action='store', type=str, dest='method', metavar='STR', default='blast', required=False,
        help='vOTU clustering method, including "blast" and "cdhit". [default=blast]'
    )

#---------------------binning------------------------
    subpsr = subparsers.add_parser(
        'binning',
        help='Cluster input viral contigs based on sequence composition, coverage and other features to generate potential viral MAGs'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr.add_argument(
        '-m', '--method', action='store', type=str, dest='method', metavar='STR', default='blast', required=False,
        help='vOTU clustering method, including "blast" and "cdhit". [default=blast]'
    )

#---------------------classify------------------------
    subpsr = subparsers.add_parser(
        'classify',
        help='Classify the viral contigs by a series of tools'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)

#---------------------compare-------------------------
    subpsr = subparsers.add_parser(
        'compare',
         help='Compare the virus protein sequence by vContact2'
    )
    subpsr = addFaArg(subpsr,'prot')
    subpsr = addGlbArg(subpsr)
    subpsr.add_argument(
        '-e', '--prefix', action='store', type=str, dest='orfprefix', metavar='STR', default=False, required=False,
        help='the first element of array which was generated by splitting ORFs ID'
    )

#---------------------vir_quant-----------------------
    subpsr = subparsers.add_parser(
        'vir_quant',
        help='Calculate the abundance and diversity of each microbial community'
    )
    #subpsr = addAllFqArg(subpsr)
    subpsr = addSampArg(subpsr)
    subpsr = addFaArg(subpsr)
    subpsr = addCheckVArg(subpsr)
    subpsr = addGlbArg(subpsr)#,1)
    subpsr = addTaxaArg(subpsr)#,1)
    subpsr.add_argument(
        '-m', '--coverm_method', action='store', type=str, dest='coverm_method', metavar='STR', default='mean', required=False,
        help='Abundance claculation method by CoverM, including "mean", "metabat" and so on. If preparing the input coverage file for VirCraft binning module, "metabat" should be chosen [default=mean]'
    )

    
#-----------------------gene_quant---------------------
    subpsr = subparsers.add_parser(
        'gene_quant',
        help='Calculate the abundance for each microbial communinity'
    )
    #subpsr = addAllFqArg(subpsr)
    subpsr = addSampArg(subpsr)
    subpsr = addFaArg(subpsr, 'gene')
    subpsr = addGlbArg(subpsr)#,1)

#-----------------------func_annot---------------------
    subpsr = subparsers.add_parser(
        'func_annot',
        help='Gene annotation and quantification'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)

#-----------------------host_pred---------------------
    subpsr = subparsers.add_parser(
        'host_pred',
        help='Predict the linkages between virus and hosts.'
    )
    subpsr = addFaArg(subpsr)
    subpsr = addGlbArg(subpsr)
    subpsr = addTaxaArg(subpsr)#,1)
    subpsr.add_argument(
        '-m', '--host_mags', action='store', type=str, dest='hostsdir', metavar='STR', default=False, required=True,
        help='A directory stored the MAGs of hosts. [required!]'
    )
    subpsr.add_argument(
        '-g', '--gtdbtk', action='store', type=str, dest='gtdbtkdir', metavar='STR', default=None,
        help='The gtdbtk results directory which performed based on host contigs. [default=None]'
    )
    return parser

#------------Functions for adding Arguments-----------
def addGlbArg(psr):
    ThreadsHelpDict = {
        0 : 'Number of processes/threads to use [default=8]',
        1 : 'Number of processes/threads to use for each task in a batch [default=8]'
    }
    psr.add_argument(
        '-d', '--config-file', action='store', type=str, dest='config', metavar='STR', default=False, required=False,
        help='Configure file can point to the parameters of certain tools and the database locations for VirCraft [default=False]'
    )
    psr.add_argument(
        '-t', '--threads', action='store', type=str, dest='threads', metavar='INT', default=8,
        help=ThreadsHelpDict[0]#batch]
    )
    psr.add_argument(
        '-u', '--unrun', action='store_true', dest='unrun', default=False, required=False,
        help='This parameter is mainly used for debugging. If this parameter is set, the script will not run directly, but will generate scripts for each analysis step [default=False]'
    )
    psr.add_argument(
        '-r', '--clear', action='store_true', dest='clear', default=False, required=False,
        help='Remove intermediate result files generated during program execution to save storage space [default=False]'
    )
    psr.add_argument(
        '-o', '--outdir', action='store', type=str, dest='outdir', metavar='STR', default=False, required=True,
        help='Output folder [default="."]'
    )
    return psr

def addProcArg(psr, dflt: str):
    HelpDict = {
        'fuc': 'Select the optional analysis process of read_qc (f, u, and/or c), i.e. "-p fuc". Among these, "f" means filter, "u" means removing the duplications and get the unique reads, and "c" refers to the process of remove the contamination from a customized reference database [default="fu"]',   
        'ms': 'Select the optional analysis process of assembly (s and/or c), i.e. "-p ms". Among these, "m" and/or "s" represent the assembly tool of MEGAHIT and/or SPAdes. i.e. "ms" refer to the process as follows: 1) assemble the reads to metagenome using MEGAHIT, 2) map all reads back to the assembled contigs and get the unmapped reads, 3) assemble the unmapped reads with SPAdes, and 4) merge the assembly results from 2) and 3) [default="ms"]'
    }
    psr.add_argument(
        '-p', '--process', action='store', type=str, dest='process', metavar='STR', default=dflt,
        help=HelpDict[dflt]
    )
    return psr

def addPairFqArg(psr):
    psr.add_argument(
        '-1', '--fastq1', action='store', type=str, dest='fq1', metavar='STR', default=False,
        required=True,help='FastQ file for read 1'
    )
    psr.add_argument(
        '-2', '--fastq2', action='store', type=str, dest='fq2', metavar='STR', default=False,
        help='FastQ file for read 2'
    )
    return psr

def addAllFqArg(psr):
    psr.add_argument(
        '-q', '--fastqs', action='store', type=str, dest='fastqs', metavar='STR', default=False,
        help='All clean FastQs, i.e. path/*.fastq'
    )
    return psr

def addFaArg(psr, seq_type='ctg'):
    HelpDict = {
        'ctg': 'The absolute or relative path to a FastA file containing viral configs or vOTUs sequences. e.g., viral_positive_contigs.fsa',
        'gene': 'The absolute or relative path to a FastA file containing gene sequences. e.g., viral_positive_contigs.ffn',
        'prot': 'The absolute or relative path to a FastA file containing protein sequences. e.g., viral_positive_contigs.faa'
    }
    psr.add_argument(
        '-a', '--fasta', action='store', type=str, dest='fasta', metavar='STR', default=False,
        help=HelpDict[seq_type]
    )
    return psr

def addSampArg(psr):
    psr.add_argument(
        '-s', '--sampinfo', action='store', type=str, dest='samp_info', metavar='STR', default=False, required=True,
        help='Sample information file with the header of \"#Sample\\tGroup\\tDataPath\\n\", and the format of each line in the text is \"sample name\\tgroup name\\tfull path of fastq1, full path of fastq2\\n\"'
    )
    return psr 

def addCutoffArg(psr):
    psr.add_argument(
        '-l', '--cutoff', action='store', type=int, dest='cutoff', metavar='INT', default=1500,
        help='The minimal length of contigs/scaffolds. [default=1500]'
    )
    return psr

def addCheckVArg(psr):
    psr.add_argument(
        '-c', '--checkv', action='store', type=str, dest='checkv', metavar='STR', default=False,
        help='CheckV results directory'
    )
    return psr

def addTaxaArg(psr):
    psr.add_argument(
        '-x', '--taxa', action='store', type=str, dest='taxa', metavar='STR', default=None,
        help='Viral taxonomic annotation file, which generated by classify module in VirCraft packages. [default=None]'
    )
    return psr
