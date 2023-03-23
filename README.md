# VirCraft
VirCraft是一个自动化病毒组学分析流程软件。
VirCraft is an automatic viromic analysis pipeline.

## 1 软件结构
![Overall workflow of VirCraft](docs/Overall_workflow_of_VirCraft.png)
#### 1.1 测序数据质控
reads_qc
#### 1.2 宏基因组组装
assembly
#### 1.3 病毒鉴定
identify
###### 1.3.1 What_the_Phage流程

[表1-1 What_the_Phage识别病毒contigs的标准参数](https://doi.org/10.1038/s42003-022-04027-y)
|tool|criteria|filter|
|:----:|:----:|:----:|
|marvel|probability according to Random Forest algorithm|>75%|
|VirFinder|p-value|>0.9|
|PPP-Meta|contig classification|"Phage"|
|VirSoter and VirSorter_virome|Category of detection (1, 2 or 3: intact, incomplete or questionable)|Category 1 & 2|
|MetaPhinder & MetaPhinder-own-DB| A) contig classification & B) average nucleotide identity % |A) Phage & B) > 50|
|DeepVirFinder|p-value|> 0.9|
|Vibrant & Vibrant_virome|contig classification|Virus|
|Phigaro|Indicator function||
|Virnet|p-value (as median across all hits per contig)|>0.5|
|Virsorter 2|dsDNA phage score|>0.9|
|Seeker|Score|>0.75|

###### 1.3.2 vir-id-sop流程

vir-id-sop流程主要依据[Guo等](dx.doi.org/10.17504/protocols.io.bwm5pc86)提供的病毒鉴定标准分析流程（Viral sequence identification SOP with VirSorter2 V.3）开发，主要包括以下步骤。
1.病毒序列鉴定。设置cutoff值为0.5，以最大的灵敏度运行VirSorter2，一般深海病毒组学项目只针对噬菌体dsDNA和ssDNA噬菌体。选择最小序列长度5000 bp。"-j"选项为CPU核数。注意，"--keep-original-seq"选项保留了环状和(接近)完整病毒contigs的原始序列(整个序列的评分为>0.8)，
2.病毒序列质控和修剪。应用checkV对VirSorter2结果进行质控，以修剪末端留下的可能的宿主基因，并处理环状contigs的重复片段。"-t"用于调整使用的CPU核数。
3.VirSorter2(>=2.2.1)再次运行。输入文件为checkv修剪的序列，结果生成的"affi-contigs.tab"文件是DRAMv识别AMG所需的文件。"--seqname-suffix-off"选项保留原始输入序列名称，因为在第二步中不可能从同一个contig中获得1个前病毒，而"--viral-gene-rich-off"选项关闭了病毒基因多于宿主基因的要求，以确保VirSorter2在这一步不进行任何筛选。
4.基因注释。DRAMv对识别的序列进行基因注释，可用于手动判定。
5.病毒序列判定。本步骤所有程序均为in-house。
```
#step 1 viral sequence identification.
virsorter run --keep-original-seq -i 5seq.fa -w vs2-pass1 --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all
#step 2 Quality Control
checkv end_to_end vs2-pass1/final-viral-combined.fa checkv -t 28 -d /fs/project/PAS1117/jiarong/db/checkv-db-v1.0
cat checkv/proviruses.fna checkv/viruses.fna > checkv/combined.fna
#step 3 Prepare for DRAMv
virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i checkv/combined.fna -w vs2-pass2 --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 28 all
#step 4 annotation
DRAM-v.py annotate -i vs2-pass2/for-dramv/final-viral-combined-for-dramv.fa -v vs2-pass2/for-dramv/viral-affi-contigs-for-dramv.tab -o dramv-annotate --skip_trnascan --threads 28 --min_contig_size 1000
DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill
#step 5 curation
sed '1s/seqname/contig_id/' vs2-pass1/final-viral-score.tsv > curation/final-viral-score.tsv #add the header for "final-viral-score.tsv" file
linkTab.py curation/final-viral-score.tsv checkv/contamination.tsv left contig_id curation/curation_vs2_checkv.tsv #combine two tables (curation/final-viral-score.tsv and checkv/contamination.tsv) by "contig_id" field.
vCurator.py . #Generate the auto-curated viral contigs and manu-curate contigs table.
cut -f 1 curation/curated_contigs.xls |grep -v "contig_id" > curation/contigs_id.list #Get the curated contigs list
sed 's/_1 / /' checkv/combined.fna > checkv/combined_modi.fna #Rename contig id
extrSeqByName.pl curation/contigs_id.list checkv/combined_modi.fna curation/virus_positive.fna #
```

#### 1.4 病毒分类
classify

#### 1.5 病毒物种丰度分析
vir_quant

#### 1.6 基因功能注释
func_annot

#### 1.7 病毒宿主分析
host_prid


## 2 软件安装和数据库部署


```
sh install.sh
```
## 3 软件使用方法
当所有依赖的软件和数据库准备就绪，VirCraft使用起来就比较简单了。VirCraft主程序脚本（virCraft.py）不包含任何功能模块，使用者可以单独使用这些模块。
```
virCraft.py -h
usage:
        ./virCraft.py -h [<options>] -o <outdir>
        subcommands: an optional functional module, including assembly, identify, votus, classify, compare, vir_quant, func_annot and host_prid.
        options: options described below in the section of Options.
        outdir: output directory.

VirCraft is an flexible pipeline for metaviromic data analysis.

optional arguments:
  -h, --help            show this help message and exit

subcommands:
  valid subcommands

  {reads_qc,assembly,identify,votus,classify,compare,vir_quant,func_annot,host_prid}
    reads_qc            Pair-end FastQ reads qualitiy control.
    assembly            Assemble the reads to contigs or scaffolds using
                        MegaHit and/or SPAdes
    identify            identify the viral contigs from a assembly fasta,
                        using vir-id-sop
    votus               construct the non-redundant virus operational
                        taxonomic unit (vOTU) reference
    classify            classify the virus contigs by Demovir
    compare             Compare the virus protein sequence by vContact2
    vir_quant           Calculate the abundance and diversity of each
                        microbial community                                        func_annot          Gene annotation and quantification
    host_prid           Predict the hosts of virus
```
各个模块可单独运行。例如，运行identify（病毒鉴定）模块。
```
usage: ./virCraft.py identify [<options>] -o <outdir>
        subcommands: an optional functional module, including assembly, identify, votus, classify, compare, vir_quant, func_annot and host_prid.
        options: options described below in the section of Options.
        outdir: output directory. identify
        [-h] [-a STR] [-t INT] -o STR

optional arguments:
  -h, --help            show this help message and exit
  -a STR, --fasta STR   The absolute/relative path of a vOTUs FastA file
  -t INT, --threads INT   Number of processes/threads to use
  -o STR, --outdir STR  Output directory
```
## 4 结果文件说明

## 5 注意事项
当前版本只能生成脚本并不能直接运行，请生成脚本后自行运行。
## 6 版本更新日志

**VirCraft-v0.0.1版**
```
初始版本。
```

**VirCraft-v0.0.2版**
```
初始版本。
```

**VirCraft-v0.0.3版**
```
尝试多路分析失败的版本。
```

**VirCraft-v0.0.4版**
```
1.暂时不支持多线程的版本。
2.quantify模块完成，包括统计多样性分析、散点图、柱状图等。
3.大势所趋，用import argparse代替from optparse import OptionParser。
```

**VirCraft-v0.0.5版**
```
添加AMGs分析，暂定用dramv和vibrant分析，尚未整合。
```

**VirCraft-v0.0.6版**
```
1.classify工具的更新：增加blast aginst NCBI viral RefSeq分类部分。
2.手动curation部分的自动化。
3.整合compare结果和classify结果。
4.完善host_pred模块。
5.部分模块支持RUN功能。
6.软件结构更新：data、general等模块。
```
