from ..config import setVari, conf
from ..process import cmdExec, general

def calcTPM(samp: str, sort_bam: str):
    envs = setVari.selectENV('VirCraft')
    tpm_cmd = [envs]
    tpm = f'{wkdir}/{samp}.tpm'
    tpm_cmd.extend(
        ['coverm contig', '-b', sort_bam, 
         '-t 20 --min-read-aligned-length 50 --min-read-percent-identity 0.95 --proper-pairs-only -m tpm', 
         '>', tpm]
    )
    tpm_sh = f'{wkdir}/{samp}_tpm.sh'
    general.printSH(bwa_idx_sh, tpm_cmd)
    results = cmdExec.execute(tpm_cmd)
    return results

def TPMBySamp(config: str, outdir: str):
    global wkdir
    groups, confDict, sampDict = conf.prepInfo(config)
    wkdir = f'{outdir}/05.abundance/2.coverm'
    bam_dir = f'{outdir}/05.abundance/1.bwa'
    general.mkdir(wkdir)
    for samp in sampDict.keys():
        sort_bam = f'{bam_dir}/{samp}.sort.bam'
        results += calcTPM(samp, sort_bam)
    merge_tpm_cmd = ['merge_tpms.pl', confDict['SampInfo'], wkdir, '\n']
    tax_anno = f'{outdir}/04.classify/DemoVir_assignments.txt'
    modi_tax_anno = f'{wkdir}/DemoVir_assignments.modi.txt'
    merged_tpm = f'{wkdir}/all_merged.tpm'
    merged_anno_tpm = f'{wkdir}/all_merged.anno.tpm'
    merged_anno_modi_tpm = f'{wkdir}/all_merged.anno.modi.tpm'
    merge_tpm_cmd.extend(
        ["sed '1s/Sequence_ID/Contig/'", tax_anno, '>', modi_tax_anno, '\n',
         'linkTab.py', merged_tpm, modi_tax_anno, 'left Contig', merged_anno_tpm, '\n',
         'tpmAddSource.py', merged_anno_tpm, merged_anno_modi_tpm, '\n',
         'pheatmap_for_tpm.R', merged_anno_modi_tpm, confDict['SampInfo'], wkdir, '\n']
    )
    quality_summary = f'{outdir}/03.vOTUs/merged/2.checkv/quality_summary.tsv'
    len_sum_tpm_qual_xls = f'{wkdir}/contig_quality_summary.xls'
    merge_tpm_cmd.extend(
         ['sumAbundance.py', merged_anno_modi_tpm, quality_summary, wkdir, '\n', 'fa_length_tpm_scatter.R', len_sum_tpm_qual_xls, wkdir, '\n']
    )
    tax_tpm = f'{wkdir}/tax_tpm.xls'
    merge_tpm_cmd.extend(
        ['abundByTax.py', merged_anno_modi_tpm, wkdir, '\n',
         'barplot_for_taxa_tpm.R', tax_tpm, wkdir, '\n']
    )
    merged_tpm_sh = f'{wkdir}/merged_anno_tpm.sh'
    general.printSH(merged_tpm_sh, merge_tpm_cmd)
    results = cmdExec.execute(merge_tpm_cmd)
    return results
