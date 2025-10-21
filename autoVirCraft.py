#!/usr/bin/env python3
#coding=utf-8
#author: yangming, yangm@idsse.ac.cn
#Created on Tus Jan 2 15:38:37 2024
import sys
from crafts.config import arguments
from crafts.config import config
from crafts.general import utils

version = '0.0.15'
parser = arguments.autOpts(sys.argv[0],version)
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()

utils.mkdir(outdir)
def generate_shell(args):
    samp_info = args.samp_info
    threads = args.threads
    host_mags = args.host
    outdir = args.outdir
    CONFIG = config.VirCfg()
    groups, sampDict = CONFIG.readSampInfo(samp_info)
    threads_per_samp = threads/len(groups)
    shell_dir = f'{outdir}/shell'
    utils.mkdir(shell_dir)
    identification_dir = f'{outdir}/01.Identification'
    vOTUs_dir = f'{outdir}/02.vOTUs'
    utils.mkdir(vOTUs_dir)
    vCTGs_dir = f'{vOTUs_dir}/viral_contigs'
    utils.mkdir(vCTGs_dir)
    ln_vCTGs_cmd = ['mkdir', vCTGs_dir]
    for samp in sampDict.keys(): #Generate the shell scripts for reads_qc
        fq1, fq2 = sampDict[samp][1].split(',')
        samp_outdir = f'{identification_dir}/{samp}'
        cmd = ['virCraft.py reads_qc', '-1', fq1, '-2', fq2, '-t', threads_per_samp, '-p fuc', '-o', samp_outdir, '\n']
        shell = f'{shell_dir}/00.{samp}_reads_qc.sh'
        utils.printSH(shell, cmd)
        clean_fq1 = f'{samp_outdir}/'
        clean_fq2 = f'{samp_outdir}/'
        cmd = ['virCraft.py assemble', '-1', clean_fq1, '-2', clean_fq2, '-t', threads_per_samp, '-p sm -l 2000', '-o', samp_outdir, '\n']
        shell = f'{shell_dir}/01.{samp}_assemble.sh'
        utils.printSH(shell, cmd)
        metagenome = f'{assemble_dir}/{sample}/final_assembly.filt.gt2000.fa'
        cmd = ['virCraft.py identify', '-a', metagenome, '-t', threads_per_samp, '-o', samp_outdir, '\n']
        shell = f'{shell_dir}/02.{samp}_identify.sh'
        utils.printSH(shell, cmd)
        ori_viral_ctg = f'{samp_outdir}/viral_positive_ctg.fna'
        link_viral_ctg = f'{vCTGs_dir}/{samp}.fa'
        link_cmd = f'&& ln {ori_viral_ctg} {link_viral_ctg}'
        ln_vCTGs_cmd.append(link_cmd)
    
    merged_vctg_fasta = f'{vOTUs_dir}/all_viral_ctg.fa'
    cmd = ln_vCTGs_cmd
    cmd.extend(['\nrenameMergeSeq.pl', '-i', vCTGs_dir, '-o', merged_vctg_fasta, '\n',
                'virCraft.py votus', '-a', ori_viral_ctg, '-o', vOTUs_dir, '\n')
    shell = f'{shell_dir}/03.votu_cluster.sh'
    utils.printSH(shell, cmd)
    votu_fasta = f'{vOTUs_outdir}/all_viral_ctg_votus.fa'
    classify_dir = f'{outdir}/04.Classification'
    cmd = ['virCraft.py classify', '-a', votu_fasta, '-t', threads, '-o', classify_dir, '\n']
    shell = f'{shell_dir}/04.votu_classify.sh'
    utils.printSH(shell, cmd)
    viral_hosts_dir = f'{outdir}/05.ViralHosts'
    host_mag_dir = 
    viral_taxon = 
    cmd = ['virCraft.py host_pred', '-a', votu_fasta, '-m', host_mag_dir, '-x', viral_taxon, '-t', threads, '-o', viral_hosts_dir, '\n']
    shell = f'{shell_dir}/05.host_predict.sh'
    utils.printSH(shell, cmd)
    viral_abundance_dir = f'{outdir}/06.ViralAbundance'
    cmd = ['virCraft.py vir_quant', '-s', samp_info, '-a', votu_fasta, '-t', threads, '-x', viral_taxon, '-c', checkv, '-m mean', '-o', viral_hosts_dir, '\n']
    shell = f'{shell_dir}/06.viral_abundance.sh'
    utils.printSH(shell, cmd)
    func_annot_dir = f'{outdir}/07.FunctionAnnotation'
    cmd = ['virCraft.py func_annot', '-a', votu_fasta, '-t', threads, '-o', func_annot_dir, '\n']
    shell = f'{shell_dir}/07.viral_gene_annotation.sh'
    utils.printSH(shell, cmd)
    gene_abundance_dir = f'{outdir}/08.GeneAbundance'
    viral_genes_ffn = f'{func_annot_dir}/prodigal/all_viral_ctg_votus.ffn'
    cmd = ['virCraft.py gene_quant', '-s', samp_info, '-a', viral_genes_ffn, '-t', threads, '-o', gene_abundance_dir, '\n']
    shell = f'{shell_dir}/08.gene_abundance.sh'
    utils.printSH(shell, cmd)
    return 0

def execute_shell():




