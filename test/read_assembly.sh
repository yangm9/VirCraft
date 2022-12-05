export "PATH=/data_backup/software/miniconda3/envs/VirCraft/bin:$PATH"
conda activate VirCraft
spades.py --pe1-1 read_1.fq --pe1-2 read_2.fq -t 8 -o /home/yangming/VirCraft/test/spades --careful -m 1300 -k 21,33,55,77,99,127 
bwa index -a bwtsw /home/yangming/VirCraft/test/spades/scaffolds.fasta -p /home/yangming/VirCraft/test/alignment/scaffoldsIDX 
bwa mem -t 8 /home/yangming/VirCraft/test/alignment/scaffoldsIDX read_1.fq read_2.fq |grep -v NM:i:> /home/yangming/VirCraft/test/alignment/unused_reads.sam 
sam_to_fastq.py /home/yangming/VirCraft/test/alignment/unused_reads.sam > /home/yangming/VirCraft/test/alignment/unused_reads.fq 
megahit -r /home/yangming/VirCraft/test/alignment/unused_reads.fq -o /home/yangming/VirCraft/test -t 8 -m 80000000000 --tmp-dir /home/yangming/VirCraft/test/megahit/megahit.tmp  cat /home/yangming/VirCraft/test/spades/scaffolds.fasta /home/yangming/VirCraft/test/megahit/final.contigs.fa > /home/yangming/VirCraft/test/final_assembly.fasta 
assemb_stat.pl /home/yangming/VirCraft/test/final_assembly.fasta /home/yangming/VirCraft/test/final_assembly.fasta >/home/yangming/VirCraft/test/stat.tab
SeqLenCutoff.pl /home/yangming/VirCraft/test/final_assembly.fasta /home/yangming/VirCraft/test/filter/scaffolds.filt 5000 
