# This dockerfile uses the centos image
# VERSION 1 - EDITION 1
# Author: zhang3
# Command format: Instruction [arguments / command] ..
FROM centos
MAINTAINER zhang3 zhang3@qq.com
RUN yum install -y wget
RUN cd /opt; \
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh; \
sh Miniconda3-latest-Linux-x86_64.sh -p /opt/miniconda3 -b
RUN /opt/miniconda3/bin/conda init
RUN source /root/.bashrc
RUN mkdir /opt/alphago/
