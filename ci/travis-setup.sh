#!/bin/bash
set -euo pipefail

curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda
export PATH=~/anaconda/bin:$PATH

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
#conda config --add channels r
conda config --add channels bioconda

conda install -y sphinx "python=3.5" "conda<4.3"

conda install -y --file requirements.txt -vv

# we also want conda-build to build a conda package
conda install -y conda-build

# for docs
pip install guzzle_sphinx_theme

~/anaconda/bin/python setup.py install
