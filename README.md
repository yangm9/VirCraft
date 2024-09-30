# VirCraft - a flexible bioinfomatic pipeline for viromic analysis from metagenomic data
VirCraft is designed to be an **easy-to-use Viral Metagenomic Analysis Craft** that integrates a variety of tools to cover all major steps in viral metagenomic analysis. Installation is straightforward by running the `setup_env` and `setup_db` modules included with VirCraft. For viral metagenomic analysis, VirCraft offers multiple modules for various tasks, including read quality control (`reads_qc` module), viral assembly (`assemble` module), identification of viral contigs (`identify` module), clustering of viral operational taxonomic units (vOTUs) (`votus` module), viral classification (`classify` module), dataset comparison (`compare` module), viral-host linkage analysis (`host_pred` module), functional annotation (`func_annot` module), and viral/gene abundance analysis (`vir_quant` and `gene_quant` modules).
VirCraft supports execution on both local and cloud-based Linux systems, utilizing containerization technology to ensure reproducibility. While there is no single best method for processing metagenomic data, VirCraft offers a fast and simple approach, ideal for initial analysis before exploring deeper parameterization. It can be applied to various environments, such as gut, water, and soil microbiomes (see the VirCraft paper for benchmarks). Each module in VirCraft is standalone, allowing users to select only the tools relevant to their data.

![Overall workflow of VirCraft](docs/Overall_workflow_of_VirCraft.png)

## VIRCRAFT MODULES

Installation modules:
1) setup_env:  conda environment installation
2) setup_db:   bioinfomatic databases deployment

#### Metagemonic data pre-processing modules

3) reads_qc:   read trimming and contaminations (e.g. human) removal
4) assemble:   metagenomic assembly and QC with metaSPAdes and/or MegaHit

#### viral sequence processing modules

5) identify:    viral contig identification with VirSorter2, VIBRANT, DeepVirFinder and CheckV
6) votus:       vOTUs clustering at species level
7) classify:    conservative but accurate taxonomic annotations for vOTUs
8) compare:     compare vOTUs with other viral datasets
9) host_pred:   viral-host linkages prediction
10) func_annot: functionally annotate genes in a set of votus or viral contigs
11) vir_quant:  estimate viral abundance and execute the diversity analysis across sample by group
12) gene_quant: estimate gene abundance

## SYSTEM REQUIREMENTS
The resource requirements for VirCraft depend heavily on the size of the data being processed. Due to the high memory demands of certain tools (like metaSPAdes), it is recommended to have at least **8 CPU cores** and **64 GB** of RAM for optimal performance. Officially, VirCraft supports **Linux x64** systems, but it can also be installed on **macOS** manually or through **Docker** for flexibility.

## INSTALLATION

#### Manual installation (recommended):
0. Install mamba: `conda install -y mamba`. Mamba will efficiently replace Conda, performing the same tasks but much faster.
1. Download or clone this ripository: `git clone https://github.com/yangm9/VirCraft.git`
2. Install the conda environments for VirCraft: `mkdir vc_envs && /your_path/VirCraft/virCraft.py setup_env -o tmp_envs`
*Note: 1) The "vc_envs" directory will store some certain tools and log files; 2) if the installation fails, simply rerun the installation commands until the all conda environments are installed; 3) Additionally, users can also manage the version of the tools in VirCraft by modifing the YAML files located in the directory: "/your_path/VirCraft/crafts/setup/conda_env_yaml/".*

Finally, a total of 12 conda environments will be installed as follows.

```
VC-Assembly              /your_path/miniconda3/envs/VC-Assembly
VC-CheckV                /your_path/miniconda3/envs/VC-CheckV
VC-DRAMv                 /your_path/miniconda3/envs/VC-DRAMv
VC-DeepVirFinder         /your_path/miniconda3/envs/VC-DeepVirFinder
VC-GTDBTk                /your_path/miniconda3/envs/VC-GTDBTk
VC-General               /your_path/miniconda3/envs/VC-General
VC-Quantify              /your_path/miniconda3/envs/VC-Quantify
VC-ReadsQC               /your_path/miniconda3/envs/VC-ReadsQC
VC-VHMatcher             /your_path/miniconda3/envs/VC-VHMatcher
VC-VIBRANT               /your_path/miniconda3/envs/VC-VIBRANT
VC-VirSorter2            /your_path/miniconda3/envs/VC-VirSorter2
VC-vContact2             /your_path/miniconda3/envs/VC-vContact2
```

*The total size of all conda environments will be **~18G**.*
(Optional) if the installation of any environment failed, users can entry it and install them manually, but the environment name should be correct.

3. Deposit the bioinformatic databases for VirCraft: `mkdir vc_db && /your_path/VirCraft/virCraft.py setup_db -o vc_db`
*Note: the "vc_db" will be an important directory to deposit the intermediate files for installing conda environment, so this directory can be deleted after the installation is complete.*
After VirCraft database installation, the size of each database will be as follows:

```
7.9G    vc_db/VC-dbCANDB
199G    vc_db/VC-iPHoPDB
112M    vc_db/VC-ContaminantsDB
30G     vc_db/VC-GTDBTkDB
166G    vc_db/VC-DRAMvDB
48G     vc_db/VC-eggNOGDB
11G     vc_db/VC-VIBRANTDB
620M    vc_db/VC-ViralRefSeqDB
11G     vc_db/VC-VirSorter2DB
6.7G    vc_db/VC-KofamScanDB
1.1G    vc_db/VC-DemovirDB
114M    vc_db/VC-DeepVirFinder
6.4G    vc_db/VC-CheckVDB
```
*The total size of all databases will be **~500G**.*

4. (Optional) In the step 2 and 3, the user can install the conda environment by manually running these scripts by adding the "-u" parameter option on the command line. If the user chooses this manual method, the config file ("/your_path/VirCraft/config") should be modified manually. Alternatively, users can also install all the conda environments by phisically transplant, and the steps are shown as follows:

1\) Download the VirCraft Conda environment packages from [Google Drive](https://drive.google.com/drive/folders/1QdwwcWSQv0IZNNFBvcENpIUWXY4adQGa?usp=sharing). For users in mainland China, download the Conda environment from Baidu Cloud: [VirCraft_conda_ENVs](https://pan.baidu.com/s/1NlRJVgYTCLYhGX1E6eE_Lw) (Access code: 6rq5)

2\) Install conda-pack: `conda install conda-pack`.

3\) Upload the downloaded conda environment packages to your server. For each conda environment, you should install it as follows:
```
mkdir ~/miniconda/envs/VC-XXXXXX  #VC-XXXXXX is the name of a certain conda environment, i.e. VC-Assembly
tar -xzf my_env.tar.gz -C ~/miniconda3/envs/VC-XXXXXX
cd ~/miniconda3/envs/VC-XXXXXX
conda-unpack # Repair environmental path
```
#### Docker installation (not recommended):
0. Prepare the bioinformatic databases for VirCraft: `mkdir vc_db && /your_path/VirCraft/virCraft.py setup_db -o vc_db`
1. docker build -t vircraft:latest -f Dockerfile .

## DETAILED PIPELINE WALKTHROUGH

![Detailed workflow of VirCraft](docs/Detailed_workflow_of_VirCraft.png)

## USAGE

#### After Manual installation
```
virCraft.py -h
usage:
        ./virCraft.py -h [<options>] -o <outdir>
        subcommands: an optional functional module, including assembly, identify, votus, classify, compare, vir_quant, func_annot and host_prid.
        options: options described below in the section of Options.
        outdir: output directory.

VirCraft is an flexible pipeline for viral metagenomic data analysis.

optional arguments:
  -h, --help            show this help message and exit

subcommands:
  valid subcommands

  {reads_qc,assembly,identify,votus,classify,compare,vir_quant,func_annot,host_prid}
    reads_qc            Pair-end FastQ reads qualitiy control
    assembly            Assemble the reads to contigs or scaffolds using MegaHit and/or SPAdes
    identify            identify the viral contigs from a assembly fasta, using vir-id-sop
    votus               construct the non-redundant virus operational taxonomic unit (vOTU) reference
    classify            classify the virus contigs by Demovir
    compare             Compare the virus protein sequence by vContact2
    vir_quant           Calculate the abundance and diversity of each microbial community
    func_annot          Gene annotation and quantification
    host_prid           Predict the hosts of virus
```

Each module is run separately. For example, to run the  module:

```
./virCraft.py identify -h
usage: /home/yangming/.vc/virCraft.py identify [-h] [-a STR] [-d STR] [-t INT] [-u] [-r] -o STR [-l INT] [-w STR]

options:
  -h, --help show this help message and exit
  -a STR, --fasta STR     The absolute or relative path to a FastA file containing viral configs or vOTUs sequences. e.g., viral_positive_contigs.fsa
  -d STR, --config-file   STR Configure file can point to the parameters of certain tools and the database locations for VirCraft [default=False]
  -t INT, --threads INT   Number of processes/threads to use [default=8]
  -u, --unrun             This parameter is mainly used for debugging. If this parameter is set, the script will not run directly, but will generate scripts for each analysis step [default=False]
  -r, --clear             Remove intermediate result files generated during program execution to save storage space [default=False]
  -o STR, --outdir STR    Output folder [default is the current folder]
  -l INT, --cutoff INT    The minimal length of contigs/scaffolds. [default=1500]
  -w STR, --sop STR       The sop/pipeline for viral contigs identification, including "viral-id-sop" and "vs2-vb-dvf". [default=vs2-vb-dvf]

#### After Manual installation

After successfully creating the docker image for VirCraft (see above, in INSTALLATION), users can run VirCraft using the command like: `docker run --rm vircraft <module_name> -c config [options] ...`. For example, to run the identify module:

```
docker run --rm vircraft identify -d /your_path/config -a /your_path/test.fasta -t 32 -o viral_identification 
```

## Acknowledgements

Author of pipeline: [Ming YANG](yangm@idsse.ac.cn)

Institution: [Institute of Deep-sea Science and Engineering，Chinese Academy of Sciences](http://www.idsse.cas.cn/)
