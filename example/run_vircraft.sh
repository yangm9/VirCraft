../virCraft.py reads_qc -1 read_1.fq -2 read_2.fq -t 32 -o results 
../virCraft.py assembly -1 read_1.fq -2 read_2.fq -t 32 -o results
../virCraft.py identify -a scaffolds.fa -t 32 -o results 
../virCraft.py votus -a scaffolds.fa -t 32 -o results
../virCraft.py classify -a scaffolds.fa -t 32 -o results
../virCraft.py quantify -a scaffolds.fa -s sample_info.xls -t 32 -o results 
