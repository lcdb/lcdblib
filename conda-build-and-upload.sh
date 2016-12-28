#!/bin/bash
export LCDBLIB_VERSION=$(./_version.py version)
export LCDBLIB_BUILD=$(./_version.py build)
conda build conda-recipe
