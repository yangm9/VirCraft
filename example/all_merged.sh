virCraft.py reads_qc -1 SRR5707419_1.fq.gz -2 SRR5707419_2.fq.gz -t 8 -o reads_qc_dir -p fuc
virCraft.py assembly -1 SRR5707419_1.fq.gz -2 SRR5707419_2.fq.gz -p ms -l 10000 -o assembly_dir
virCraft.py identify -a input.fasta -t 8 -l 10000 -w vs2-vb-dvf -o vs2vbdvf_identify_dir
virCraft.py votus -a input.fasta -t 8 -o votus_dir
virCraft.py classify -a input.fasta -t 8 -o func_annot_dir
virCraft.py compare -a input.fasta -t 8 -o compare_dir
virCraft.py func_annot -a input.fasta -t 8 -o func_annot_dir
virCraft.py host_prid -a input.fasta -t 8 -o host_prid_dir -m gtdbtk_dir
