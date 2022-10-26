export PATH="/data_backup/software/miniconda3/bin/conda:$PATH"
conda activate VirCraft
spades.py --pe1-1 /home/yangming/VirCraft/examples/result/01.assembly/0.fastq/cariaco_1.fq --pe1-2 /home/yangming/VirCraft/examples/result/01.assembly/0.fastq/cariaco_2.fq --careful -t 30 -m 1300 -k 21,33,55,77,99,127 -o /home/yangming/VirCraft/examples/result/01.assembly/cariaco