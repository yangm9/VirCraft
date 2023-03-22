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

###### 1.3.2 VirSorter2流程

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
