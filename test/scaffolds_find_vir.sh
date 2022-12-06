export "PATH=/data_backup/software/miniconda3/bin:$PATH"
conda activate VirCraft
virsorter run --keep-original-seq -i /home/yangming/VirCraft/scaffolds.fa -d /data_backup/database/virsorter2DB/db-vs2 -w /home/yangming/VirCraft/test/vs2-pass1 --include-groups dsDNAphage,ssDNA -j 8 --min-length 5000 --min-score 0.5 all
checkv end_to_end /home/yangming/VirCraft/test/vs2-pass1/final-viral-combined.fa /home/yangming/VirCraft/test/checkv -d /data_backup/database/checkvDB/checkv-db-v1.2 -t 8 
cat /home/yangming/VirCraft/test/checkv/proviruses.fna /home/yangming/VirCraft/test/checkv/viruses.fna > /home/yangming/VirCraft/test/checkv/combined.fna 
virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i /home/yangming/VirCraft/test/checkv/combined.fna -d /data_backup/database/virsorter2DB/db-vs2 -w /home/yangming/VirCraft/test/vs2-pass2 --include-groups dsDNAphage,ssDNA -j 8 --min-length 5000 --min-score 0.5 all
DRAM-v.py annotate -i /home/yangming/VirCraft/test/vs2-pass2/for-dramv/final-viral-combined-for-dramv.fa -v /home/yangming/VirCraft/test/vs2-pass2/for-dramv/viral-affi-contigs-for-dramv.tab -o /home/yangming/VirCraft/test/dramv-annotate --threads 8 --skip_trnascan --min_contig_size 1000
DRAM-v.py distill -i /home/yangming/VirCraft/test/dramv-annotate/annotations.tsv -o /home/yangming/VirCraft/test/dramv-distill 
sed '1s/seqname/contig_id/' /home/yangming/VirCraft/test/vs2-pass1/final-viral-score.tsv > /home/yangming/VirCraft/test/curation/final-viral-score.tsv 
linkTab.py /home/yangming/VirCraft/test/curation/final-viral-score.tsv /home/yangming/VirCraft/test/checkv/contamination.tsv left contig_id /home/yangming/VirCraft/test/curation/curation_vs2_checkv.tsv 
vCurator.py /home/yangming/VirCraft/test 
cut -f 1 /home/yangming/VirCraft/test/curation/curated_contigs.xls |grep -v "contig_id" > /home/yangming/VirCraft/test/curation/contigs_id.list 
sed 's/_1 / /' /home/yangming/VirCraft/test/checkv/combined.fna > /home/yangming/VirCraft/test/checkv/combined_modi.fna 
extrSeqByName.pl /home/yangming/VirCraft/test/curation/contigs_id.list /home/yangming/VirCraft/test/checkv/combined_modi.fna /home/yangming/VirCraft/test/curation/virus_positive.fna 
