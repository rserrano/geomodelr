.. geomodelr documentation master file, created by
   sphinx-quickstart on Thu Dec  1 18:22:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geomodelr's query tool documentation!
================================================

Geomodelr is a web tool for creating geological models easily.
To create a geological model, go to https://www.geomodelr.com.
After creating your geological model you might want to use it
for calculations, geostatistics, simulations, or simply to know
what geological unit is present at a given path. With this tool 
you can do all that.

To use geomodelr query tool you just need to:

    import geomodelr
    # load your model.
    model = geomodelr.model_from_file('/path/to/your/model_version.json')
    # query your model.
    unit, distance = model.closest((1000, 1000, 0.0))
    # do stuff...
    if unit == 'Batholith':
    ...

You can also use this tool as a script.

    $ geomodelr -q /path/to/your/model_version.json
    x y z
    Batholith


Features
--------
- Query the model in the coordinate system you defined.
- Query the topography heights and query the model with topography in mind.
- Query the intersection of faults and planes.
- Generate grids, use it as a help tool to generate meshes, assign properties for simulations or create block models.
- What do you want to do?

Installation
------------
Install project by calling:

    pip install geomodelr

It needs boost libraries and C++ compiler. In case boost libraries are not in a 
standard location, call it with INCLUDE_DIRS=... and/or LIBRARY_DIRS=...

Support
-------
If you are having problems, write to support@geomodelr.com.

License
-------
This project is licensed under the Affero GPL license https://www.gnu.org/licenses/

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

