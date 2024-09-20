# This dockerfile uses the centos image
# VERSION 1 - EDITION 1
# Author: yangm9

FROM ubuntu:20.04

LABEL maintainer="yangm9 yangm9@qq.com"

ENV VOLUME_SIZE=20G #116M 12st
WORKDIR /opt

COPY . .

RUN apt update && \ #158M 120s
    # 安装git, wget
    apt install -y --no-install-recommends git wget build-essential cmake && \
    # 安装conda和mamba
    wget -qO /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh /tmp/miniconda.sh -p /opt/miniconda3 -b && \
    rm -f /tmp/miniconda.sh && \
    # 初始化conda并安装mamba
    /opt/miniconda3/bin/conda init && \
    /opt/miniconda3/bin/conda install -y -c conda-forge mamba && \ #1.3G 3m4.330s
    # 安装VirCraft
    export PATH="/opt/VirCraft:$PATH" && \ 
    virCraft.py setup_env -w -o envs && \
    conda clean --all --yes && \
    rm -rf /var/lib/apt/lists/* && \
    # 创建非root用户
    groupadd -r appuser && useradd -r -g appuser appuser

# 切换到非root用户
USER appuser

# 暴露端口
EXPOSE 518
EXPOSE 443

# 运行应用程序
ENTRYPOINT ["conda", "run", "-n", "base", "virCraft.py"]
