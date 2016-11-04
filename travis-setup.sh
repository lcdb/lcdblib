#!/bin/bash
set -euo pipefail

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
for channel in $(grep -v "^#" bioconda_utils/channel_order.txt); do
    conda config --add channels $channel
done

conda install -y --file requirements_dev.txt
python setup.py install
