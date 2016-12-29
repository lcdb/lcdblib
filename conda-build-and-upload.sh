#!/bin/bash

# Always build with conda to test any build issues, but only upload if we're on
# master branch and it's not a pull request.
export LCDBLIB_VERSION=$(./_version.py version)
export LCDBLIB_BUILD=$(./_version.py build)
conda build conda-recipe

if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
  conda install anaconda-client -y
  anaconda \
    -t $ANACONDA_TOKEN \
    upload \
    -u daler \
    /home/travis/anaconda/conda-bld/linux-64/lcdblib-${LCDBLIB_VERSION}-${LCDBLIB_BUILD}.tar.bz2
fi
