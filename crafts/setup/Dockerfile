# This dockerfile uses the centos image
# VERSION 1 - EDITION 1
# Author: zhang3
# Command format: Instruction [arguments / command] ..
FROM centos
MAINTAINER yangming yangm@idsse.ac.cn
RUN yum install -y wget && \
cd /opt && \
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
sh Miniconda3-latest-Linux-x86_64.sh -p /opt/miniconda3 -b && \
/opt/miniconda3/bin/conda init && \
source /root/.bashrc && \
mkdir /opt/VirCraft/ && \

