export "PATH=/data_backup/software/miniconda3/bin:$PATH"
conda activate VirCraft
cd-hit-est -i /home/yangming/VirCraft/scaffolds.fa -o /home/yangming/VirCraft/test/scaffolds_votus.fa -T 8 -c 0.95 -aS 0.85 -n 10 -d 0 -M 160000
checkv end_to_end /home/yangming/VirCraft/test/scaffolds_votus.fa /home/yangming/VirCraft/test/checkv -d /data_backup/database/checkvDB/checkv-db-v1.2 -t 8 
cat /home/yangming/VirCraft/test/checkv/proviruses.fna /home/yangming/VirCraft/test/checkv/viruses.fna > /home/yangming/VirCraft/test/checkv/combined.fna 
