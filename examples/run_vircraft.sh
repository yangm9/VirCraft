../virCraft.py reads_qc 
../virCraft.py assembly -1 test1.fq -2 test2.fq -t 32 -o results
../virCraft.py identify -c vircraft.cfg -o results
../virCraft.py votus -c vircraft.cfg -o results
../virCraft.py classify -c vircraft.cfg -o results
