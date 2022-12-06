export "PATH=/data_backup/software/miniconda3/bin:$PATH"
conda activate VirCraft
ln -s /data_backup/software/biotools/demovir/TrEMBL_viral_taxa.RDS /home/yangming/VirCraft/test/TrEMBL_viral_taxa.RDS 
ln -s /data_backup/software/biotools/demovir/uniprot_trembl.viral.udb /home/yangming/VirCraft/test/uniprot_trembl.viral.udb 
cp /home/yangming/VirCraft/bin/demovir.* /home/yangming/VirCraft/test 
/home/yangming/VirCraft/test/demovir.sh /home/yangming/VirCraft/test/03.vOTUs/merged_virus_positive_nodup.fa 8