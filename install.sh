wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b
conda init
source ~/.bashrc
pip install -r requirements.txt
conda create -y -n VirCraft megahit -c bioconda
conda activate VirCraft
conda install -y -c bioconda trimmomatic prodigal spades fastp sortmerna bowtie2 virsorter=2
conda install -y checkv
conda install -y dram
conda install -y -c bioconda cd-hit
conda install -y -c bioconda eggnog-mapper
conda install -y -c bioconda vibrant==1.2.1
conda install -y iqtree
conda install -y -c bioconda diamond
conda install -y coverm
conda install -c bioconda crass
conda install -c bioconda -y minced
conda install -c bioconda trnascan-se
conda install -y salmon
conda install -y r-base
conda install -c bioconda hmmer
conda install -c bioconda muscle
conda install -y drep
conda install -c bioconda aragorn
