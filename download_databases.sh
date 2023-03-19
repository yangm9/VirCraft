#!/bin/bash
. "/data_backup/software/miniconda3/etc/profile.d/conda.sh"
conda activate VirCraft
virsorter setup -d db-vs2 -j 4
checkv download_database .
DRAM-setup.py prepare_databases --skip_uniref --output_dir db-dramv
