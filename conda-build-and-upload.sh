#!/bin/bash

export LCDBLIB_VERSION=$(python -c 'from lcdblib import __conda_version__; print(__conda_version__)')
export LCDBLIB_BUILD=$(python -c 'from lcdblib import __conda_build__; print(__conda_build__)')

conda build conda-recipe

# Always build with conda to test any build issues, but only upload if we're on
# master branch and it's not a pull request.
if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
  conda install anaconda-client -y
  anaconda \
    -t $ANACONDA_TOKEN \
    upload \
    -u lcdb \
    $(conda build --output conda-recipe)
fi
