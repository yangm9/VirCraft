#!/bin/bash
if conda --version >/dev/null 2>&1; then
    echo "Conda has been installed."
else
    echo "Now conda will be installed..."
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh -b
    conda init
    source ~/.bashrc
    echo "The installation of conda is complete！"
fi
#reads_qc
envs=(reads_qc assembly deepvirfinder viral-id-sop vibrant gtdbtk vContact VirCraft kofamscan)

for env in "${envs[@]}"
do
    if conda env list | grep -q $env; then
        echo "reads_qc exists!"
    else
        conda env create -f $env.yaml
        if 
    fi
done

if grep -q "DeepVirFinder" ../crafts/config/vircraft.cfg; then
    line=grep -q "DeepVirFinder" ../crafts/config/vircraft.cfg

else
    git clone https://github.com/jessieren/DeepVirFinder
    dvf_path="$(pwd)/DeepVirFinder/dvf.py"
    echo "DeepVirFinder = ${dvf_path}" >> ../crafts/config/vircraft.cfg
fi

