#!/bin/bash
set -euo pipefail
set -x

# Sets up travis-ci environment for testing bioconda-utils.

if [[ $TRAVIS_OS_NAME = "linux" ]]
then
    tag=Linux
else
    tag=MacOSX
fi
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-${tag}-x86_64.sh
sudo bash Miniconda3-latest-${tag}-x86_64.sh -b -p /anaconda
sudo chown -R $USER /anaconda
export PATH=/anaconda/bin:$PATH

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

conda install -y --file requirements.txt
