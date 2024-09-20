# This dockerfile uses the centos image
# VERSION 1 - EDITION 1
# Author: yangm9

FROM continuumio/miniconda3:latest

LABEL maintainer="yangm9 yangm9@qq.com"

ENV VOLUME_SIZE=20G 
WORKDIR /opt/VirCraft

COPY . .

RUN apt update && \ 
    # 安装git, wget
    apt install -y --no-install-recommends build-essential &&\ 
    # 初始化conda并安装mamba
    conda install -y -c conda-forge mamba && \ 
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
