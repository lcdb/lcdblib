#!/bin/bash

$PYTHON setup.py clean install --single-version-externally-managed --record /tmp/$PKG_NAME.log

# copy scripts over
mv ./bin/* "$PREFIX/bin/"
