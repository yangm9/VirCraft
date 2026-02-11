#!/usr/bin/env python3
# coding=utf-8
# author: yangming, yangm@idsse.ac.cn
# Updated by GPT-5.2-Codex

import os
import shlex
import sys
from crafts.config import arguments
from crafts.config import config
from crafts.general import utils

version = '0.0.16'


def _q(value):
    return shlex.quote(str(value))


def _step_enabled(steps: str, step_no: str) -> bool:
    return step_no in set(str(steps))


def _parse_fastq_pair(data_path: str):
    items = [x.strip() for x in data_path.split(',') if x.strip()]
    if len(items) != 2:
        raise ValueError(f'DataPath must contain exactly two FASTQ paths separated by comma: {data_path}')
    return items[0], items[1]


def generate_shell(args):
    samp_info = args.samp_info
    outdir = os.path.abspath(args.outdir)
    threads = int(args.threads)
    steps = str(args.steps)
    process = getattr(args, 'process', 'fuc')
    min_len = int(getattr(args, 'min_len', 2000))
    identify_methods = getattr(args, 'methods', 'gn')
    identify_filter = getattr(args, 'filt_mode', 'permissive')
    votu_method = getattr(args, 'votu_method', 'blast')
    coverm_method = getattr(args, 'coverm_method', 'mean')

    host_mags = getattr(args, 'hostsdir', None)
    gtdbtkdir = getattr(args, 'gtdbtkdir', None)
    taxa_input = getattr(args, 'taxa', None)
    checkv_dir = getattr(args, 'checkv', None)

    utils.mkdir(outdir)
    shell_dir = f'{outdir}/shell'
    utils.mkdir(shell_dir)

    identification_dir = f'{outdir}/01.Identification'
    vOTUs_dir = f'{outdir}/02.vOTUs'
    vCTGs_dir = f'{vOTUs_dir}/viral_contigs'
    classify_dir = f'{outdir}/04.Classification'
    viral_hosts_dir = f'{outdir}/05.ViralHosts'
    viral_abundance_dir = f'{outdir}/06.ViralAbundance'
    func_annot_dir = f'{outdir}/07.FunctionAnnotation'
    gene_abundance_dir = f'{outdir}/08.GeneAbundance'

    for d in [identification_dir, vOTUs_dir, vCTGs_dir]:
        utils.mkdir(d)

    CONFIG = config.VirCfg()
    groups, samp_dict = CONFIG.readSampInfo(samp_info)
    samp_n = max(1, len(samp_dict))
    threads_per_samp = max(1, threads // samp_n)

    sample_shs = []
    link_cmds = []

    for samp in sorted(samp_dict.keys()):
        fq1, fq2 = _parse_fastq_pair(samp_dict[samp][1])
        samp_outdir = f'{identification_dir}/{samp}'
        utils.mkdir(samp_outdir)

        cmd_lines = []

        # 1) reads_qc
        clean_fq1 = os.path.join(samp_outdir, os.path.basename(fq1).replace('.gz', ''))
        clean_fq2 = os.path.join(samp_outdir, os.path.basename(fq2).replace('.gz', ''))
        if _step_enabled(steps, '1'):
            cmd_lines.append(
                f"virCraft.py reads_qc -1 {_q(fq1)} -2 {_q(fq2)} -t {threads_per_samp} -p {_q(process)} -o {_q(samp_outdir)}"
            )
        else:
            clean_fq1, clean_fq2 = fq1, fq2

        # 2) assemble + identify
        metagenome = f'{samp_outdir}/final_assembly.gt{min_len}.fa'
        if _step_enabled(steps, '2'):
            cmd_lines.append(
                f"virCraft.py assemble -1 {_q(clean_fq1)} -2 {_q(clean_fq2)} -t {threads_per_samp} -p ms -l {min_len} -o {_q(samp_outdir)}"
            )
            cmd_lines.append(
                f"virCraft.py identify -a {_q(metagenome)} -t {threads_per_samp} -l {min_len} -m {_q(identify_methods)} -f {_q(identify_filter)} -o {_q(samp_outdir)}"
            )

        ori_viral_ctg = f'{samp_outdir}/viral_positive_ctg.fna'
        link_viral_ctg = f'{vCTGs_dir}/{samp}.fa'
        if _step_enabled(steps, '2'):
            link_cmds.append(f'ln -sf {_q(ori_viral_ctg)} {_q(link_viral_ctg)}')

        shell_path = f'{shell_dir}/00_02.{samp}.sh'
        if cmd_lines:
            utils.printSH(shell_path, ['\n'.join(cmd_lines) + '\n'])
            sample_shs.append(shell_path)

    merged_vctg_fasta = f'{vOTUs_dir}/all_viral_ctg.fa'
    if args.fasta:
        merged_vctg_fasta = os.path.abspath(args.fasta)

    votu_shell = f'{shell_dir}/03.votu_cluster.sh'
    if _step_enabled(steps, '3'):
        if _step_enabled(steps, '2'):
            cmd = [
                f'mkdir -p {_q(vCTGs_dir)}\n',
                '\n'.join(link_cmds) + '\n' if link_cmds else '',
                f'renameMergeSeq.pl -i {_q(vCTGs_dir)} -o {_q(merged_vctg_fasta)}\n',
                f'virCraft.py votus -a {_q(merged_vctg_fasta)} -t {threads} -m {_q(votu_method)} -l {min_len} -o {_q(vOTUs_dir)}\n',
            ]
        else:
            cmd = [f'virCraft.py votus -a {_q(merged_vctg_fasta)} -t {threads} -m {_q(votu_method)} -l {min_len} -o {_q(vOTUs_dir)}\n']
        utils.printSH(votu_shell, cmd)

    votu_fasta = f'{vOTUs_dir}/{os.path.splitext(os.path.basename(merged_vctg_fasta))[0]}_votus.fa'

    classify_shell = f'{shell_dir}/04.votu_classify.sh'
    if _step_enabled(steps, '4'):
        cmd = [f'virCraft.py classify -a {_q(votu_fasta)} -t {threads} -o {_q(classify_dir)}\n']
        utils.printSH(classify_shell, cmd)

    viral_taxon = taxa_input or f'{classify_dir}/{os.path.splitext(os.path.basename(votu_fasta))[0]}.votu.taxa.txt'

    host_shell = f'{shell_dir}/05.host_predict.sh'
    if _step_enabled(steps, '5') and host_mags:
        cmd = [
            f'virCraft.py host_pred -a {_q(votu_fasta)} -m {_q(host_mags)} '
            f'-x {_q(viral_taxon)} -t {threads} -o {_q(viral_hosts_dir)}'
        ]
        if gtdbtkdir:
            cmd[0] += f' -g {_q(gtdbtkdir)}'
        cmd[0] += '\n'
        utils.printSH(host_shell, cmd)

    vir_ab_shell = f'{shell_dir}/06.viral_abundance.sh'
    if _step_enabled(steps, '6'):
        cmd = [
            f'virCraft.py vir_quant -s {_q(samp_info)} -a {_q(votu_fasta)} '
            f'-t {threads} -x {_q(viral_taxon)} -m {_q(coverm_method)} -o {_q(viral_abundance_dir)}'
        ]
        if checkv_dir:
            cmd[0] += f' -c {_q(checkv_dir)}'
        cmd[0] += '\n'
        utils.printSH(vir_ab_shell, cmd)

    func_shell = f'{shell_dir}/07.viral_gene_annotation.sh'
    if _step_enabled(steps, '7'):
        cmd = [f'virCraft.py func_annot -a {_q(votu_fasta)} -t {threads} -o {_q(func_annot_dir)}\n']
        utils.printSH(func_shell, cmd)

    gene_shell = f'{shell_dir}/08.gene_abundance.sh'
    if _step_enabled(steps, '8'):
        viral_genes_ffn = f'{func_annot_dir}/prodigal/{os.path.splitext(os.path.basename(votu_fasta))[0]}.ffn'
        cmd = [f'virCraft.py gene_quant -s {_q(samp_info)} -a {_q(viral_genes_ffn)} -t {threads} -o {_q(gene_abundance_dir)}\n']
        utils.printSH(gene_shell, cmd)

    # master one-click shell
    step_shells = []
    step_shells.extend(sample_shs)
    for item in [votu_shell, classify_shell, host_shell, vir_ab_shell, func_shell, gene_shell]:
        if os.path.exists(item):
            step_shells.append(item)

    master_shell = f'{shell_dir}/run_all.sh'
    master_cmd = [f'bash {_q(sh)}\n' for sh in step_shells]
    if master_cmd:
        utils.printSH(master_shell, master_cmd)

    return step_shells, master_shell


def execute_shell(master_shell: str, unrun=False):
    if not os.path.exists(master_shell):
        print(f'No executable pipeline shell generated: {master_shell}')
        return 1
    if unrun:
        print(f'Pipeline shell generated: {master_shell}')
        return 0
    return utils.execute(master_shell)


def setup_all(args):
    outdir = os.path.abspath(args.outdir)
    threads = int(args.threads)
    unrun = args.unrun
    cmds = [
        f'virCraft.py setup_env -t {threads} -o {_q(outdir)}/envs' + (' -u' if unrun else ''),
        f'virCraft.py setup_db -t {threads} -o {_q(outdir)}/db' + (' -u' if unrun else ''),
    ]
    shell_dir = f'{outdir}/shell'
    utils.mkdir(shell_dir)
    shell = f'{shell_dir}/setup_all.sh'
    utils.printSH(shell, ['\n'.join(cmds) + '\n'])
    return execute_shell(shell, unrun=unrun)


def clear_results(args):
    outdir = os.path.abspath(args.outdir)
    rm_targets = [
        f'{outdir}/01.Identification',
        f'{outdir}/02.vOTUs/work_files',
        f'{outdir}/04.Classification/work_files',
        f'{outdir}/05.ViralHosts/work_files',
        f'{outdir}/06.ViralAbundance/work_files',
        f'{outdir}/07.FunctionAnnotation/work_files',
        f'{outdir}/08.GeneAbundance/work_files',
    ]
    for target in rm_targets:
        if os.path.exists(target):
            cmd = ['rm', '-rf', target]
            utils.show_cmd(cmd)
            os.system(' '.join(cmd))
    return 0


def main():
    parser = arguments.autOpts(sys.argv[0], version)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    if sys.argv[1] == 'exec':
        __, master_shell = generate_shell(args)
        return execute_shell(master_shell, unrun=args.unrun)
    if sys.argv[1] == 'setup':
        return setup_all(args)
    if sys.argv[1] == 'clear':
        return clear_results(args)

    parser.print_help()
    return 1


if __name__ == '__main__':
    sys.exit(main())
