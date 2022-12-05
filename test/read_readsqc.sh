export "PATH=/data_backup/software/miniconda3/envs/VirCraft/bin:$PATH"
conda activate VirCraft
ls read_1.fq read_2.fq > /home/yangming/VirCraft/test/fastp/read_list.txt 
fastuniq -i /home/yangming/VirCraft/test/fastp/read_list.txt -o /home/yangming/VirCraft/test/fastuniq/read_1.fq -p /home/yangming/VirCraft/test/fastuniq/read_2.fq 
ln -s /home/yangming/VirCraft/test/fastuniq/read_1.fq /home/yangming/VirCraft/test/read_1.fq
ln -s /home/yangming/VirCraft/test/fastuniq/read_2.fq /home/yangming/VirCraft/test/read_2.fq
