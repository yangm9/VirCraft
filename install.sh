#!/bin/bash
if conda --version >/dev/null 2>&1; then
    echo "Conda has been installed."
else
    echo "Now conda will be installed..."
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh -b
    conda init
    source ~/.bashrc
    echo "The installation of conda is complete！"
fi
conda create -y -n VirCraft megahit -c bioconda
conda activate VirCraft
conda install -y -c bioconda bwa=0.7.17 blast=2.13.0 trimmomatic=0.39 \
prodigal=2.6.3 spades=3.15.4 fastp=0.23.2 sortmerna=4.3.4 bowtie2=2.4.4 \
virsorter=2.2.3 checkv=0.8.1 dram=1.3.5 cd-hit=4.8.1 eggnog-mapper=2.1.7 \
vibrant=1.2.1 iqtree=2.2.0.3 diamond=2.0.15 coverm=0.6.1 crass=1.0.1 \
minced=0.4.2 trnascan-se=2.0.11 salmon=0.14.2 hmmer=3.3.2 muscle=5.1 \
drep=3.3.1 aragorn=1.2.41 kofamscan=1.3.0 fastuniq=1.1
conda install -c conda-forge r-base=4.2.1 r-ggplot2=3.4.0 r-vctrs=0.5.1
conda install -c conda-forge scikit-learn-intelex
pip install -r requirements.txt
pip install Bio
#pip install "scikit_learn==0.22.2.post1" #FutureWarning: The sklearn.neural_network.multilayer_perceptron module is  deprecated in version 0.22 and will be removed in version 0.24.
#pip install --upgrade scikit-learn==0.21.3
pip install --upgrade scikit-learn
## Install usearch
cd ~/VirCraft/bin
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip usearch11.0.667_i86linux32.gz
ln -s usearch11.0.667_i86linux32 usearch
## Install Demovir
git clone https://github.com/feargalr/Demovir.git

