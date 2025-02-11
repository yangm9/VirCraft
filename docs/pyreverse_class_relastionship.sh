#!/bin/bash

# 检查 pylint 是否已经安装
if ! mamba list pylint &>/dev/null; then
    echo "Install pylint..."
    mamba install -y -c conda-forge pylint
else
    echo "pylint has been installed!"
fi

# 检查 graphviz 是否已经安装
if ! mamba list graphviz &>/dev/null; then
    echo "Install graphviz..."
    mamba install -y -c conda-forge graphviz
else
    echo "graphviz has been installed!"
fi

pyreverse -o pdf -d . -p VirCraft --colorized --max-color-depth 8 ../crafts
