export "PATH=/data_backup/software/miniconda3/envs/VirCraft/bin:$PATH"
conda activate VirCraft
fastp -i read_1.fq -o /home/yangming/VirCraft/test/fastp/read_1.fq -I read_2.fq -O /home/yangming/VirCraft/test/fastp/read_2.fq -w 16 -q 20 -u 20 -g -c -W 5 -3 -l 50 -h /home/yangming/VirCraft/test/fastp/read_report.html 
ls /home/yangming/VirCraft/test/fastp/read_1.fq /home/yangming/VirCraft/test/fastp/read_2.fq > /home/yangming/VirCraft/test/fastp/read_list.txt 
fastuniq -i /home/yangming/VirCraft/test/fastp/read_list.txt -o /home/yangming/VirCraft/test/fastuniq/read_1.fq -p /home/yangming/VirCraft/test/fastuniq/read_2.fq 
bowtie2 -p 40 -N 1 -x /data_backup/database/contaminations/contaminations -l /home/yangming/VirCraft/test/fastuniq/read_1.fq -2 /home/yangming/VirCraft/test/fastuniq/read_2.fq --un-conc /home/yangming/VirCraft/test/decontaminate/read 
ln -s /home/yangming/VirCraft/test/decontaminate/read_1.fq /home/yangming/VirCraft/test/read_1.fq
ln -s /home/yangming/VirCraft/test/decontaminate/read_2.fq /home/yangming/VirCraft/test/read_2.fq
