# This dockerfile uses the centos image
# VERSION 1 - EDITION 1
# Author: yangm9
# Command format: Instruction [arguments / command] ..
FROM ubuntu
MAINTAINER yangm9 yangm9@qq.com
ENV VOLUME_SIZE=20G #116M 12st
RUN apt update && \ #158M 120s
# 安装git, wget
apt install -y git && \ #237M 30s
apt install -y wget && \ #237M 8s
apt install -y build-essential && \
apt install -y cmake && \ #237M 8s
# 安装conda和mamba
cd /opt && \ 
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \ #352M 6s
sh Miniconda3-latest-Linux-x86_64.sh -p /opt/miniconda3 -b && \ #939M 20s
/opt/miniconda3/bin/conda init && \ #939M 1s
source /root/.bashrc && \
/opt/miniconda3/bin/conda install -y -c conda-forge mamba && \ #1.3G 3m4.330s
# 安装VirCraft
git clone https://gitee.com/brightyoung/VirCraft.git && \ #1.4G 2s
export PATH="/opt/VirCraft:/opt/VirCraft/bin":$PATH && \ 
cd /opt/VirCraft && \
virCraft.py setup_env -w -o envs
#

EXPOSE 315
EXPOSE 443

#总共大概19G
