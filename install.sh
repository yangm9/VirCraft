wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b
conda init
source ~/.bashrc
pip install -r requirements.txt
conda create -y -n VirCraft megahit -c bioconda
conda activate VirCraft
conda install -y -c bioconda trimmomatic prodigal spades fastp sortmerna bowtie2 virsorter=2 checkv dram cd-hit eggnog-mapper vibrant==1.2.1 iqtree diamond coverm crass minced trnascan-se salmon hmmer muscle drep aragorn kofamscan
conda install -c conda-forge r-base
## Install usearch
cd ~/VirCraft/bin
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip usearch11.0.667_i86linux32.gz
ln -s usearch11.0.667_i86linux32 usearch
## Install Demovir
git clone https://github.com/feargalr/Demovir.git

