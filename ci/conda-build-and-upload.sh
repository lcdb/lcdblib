#!/bin/bash
set -eo pipefail

export LCDBLIB_VERSION=$(python -c 'from lcdblib import __conda_version__; print(__conda_version__)')
export LCDBLIB_BUILD=$(python -c 'from lcdblib import __conda_build__; print(__conda_build__)')

if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
  conda install anaconda-client -y
fi

# Build packages for all supported versions of Python. If we're on the master
# branch and on Travis-CI, then also upload the built package.
for PY in 35 36 37; do
    CONDA_PY=$PY conda build conda-recipe
    if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
      anaconda \
        -t $ANACONDA_TOKEN \
        upload \
        -u lcdb \
        $(CONDA_PY=$PY conda build --output conda-recipe)
    fi
done
