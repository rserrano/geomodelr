#!/bin/bash

# Build package.
# pip install -e .
# Bash should exit at any error.
set -e
# Build documentation.
cd docs
rm -rf build/rst

# RTD does not support rst builder, so have to do this.
cp source/conf.py.orig source/conf.py

# RTD can't install geomodelr, so have to do use rst builder.
cp source/geomodelr.rst.orig source/geomodelr.rst
pwd
make rst
cp build/rst/geomodelr.rst source/geomodelr.rst

# Put again without rst builder.
cp source/conf.py.rtd source/conf.py
cd ..

# Distribute linux package.
python setup.py sdist
python setup.py sdist upload

git push
