#!/bin/bash

# Build package.
pip install -e .

# Build documentation.
cd docs
touch source/geomodelr.rst
make html
cd ..

# Distribute linux package.
python setup.py sdist
python setup.py sdist upload

git push
