#!/bin/bash

# Build package.
# pip install -e .
# Bash should exit at any error.
set -e
# Build documentation.
cd docs
rm -rf build/rst
cp source/geomodelr.rst.orig source/geomodelr.rst
make rst
cp build/rst/geomodelr.rst source/geomodelr.rst
cd ..

# Distribute linux package.
python setup.py sdist
python setup.py sdist upload

git push
