# Optimized Dockerfile using Ubuntu
# VERSION 2 - EDITION 1
# Author: yangm9
# Command format: Instruction [arguments / command] ..

FROM ubuntu:20.04
MAINTAINER yangm9 yangm9@qq.com

# Set environment variables
ENV VOLUME_SIZE=20G \
    DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/miniconda3 \
    PATH=/opt/VirCraft:/opt/VirCraft/bin:$PATH

# Install system dependencies, minimize layers
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    build-essential \
    cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniconda and Mamba in one layer
RUN cd /opt && \
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh Miniconda3-latest-Linux-x86_64.sh -p $CONDA_DIR -b && \
    $CONDA_DIR/bin/conda init bash && \
    $CONDA_DIR/bin/conda install -y -c conda-forge mamba && \
    conda clean --all --yes && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Clone and install VirCraft, setup environment and database
RUN git clone https://gitee.com/brightyoung/VirCraft.git /opt/VirCraft && \
    /opt/VirCraft/virCraft.py setup_env -w -o envs && \
    /opt/VirCraft/virCraft.py setup_db -o VC-db && \
    conda clean --all --yes

# Expose required ports
EXPOSE 315 443

# Set entrypoint
ENTRYPOINT ["/bin/bash"]
