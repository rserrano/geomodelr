
Geomodelr, Introduction.
************************
 
To use geomodelr query tool you just need to::

    import geomodelr
    # load your model.
    model = geomodelr.model_from_file('/path/to/your/model_version.json')
    # query your model.
    unit, distance = model.closest((1000, 1000, 0.0))
    # do stuff...
    if unit == 'Batholith':
        ...

You can also use this tool as a script::

    $ geomodelr -q /path/to/your/model_version.json
    1000 1000 0
    Batholith


Features
========
- Query the model in the coordinate system you defined.
- Query the topography heights and query the model with topography in mind.
- Query the intersection of faults and planes.
- Generate grids, use it as a help tool to generate meshes, assign properties for simulations or create block models.
- What do you want to do?

Installation
============
Install project by calling::

    pip install geomodelr

It needs boost libraries and C++ compiler. In case boost libraries are not in a 
standard location, call it with :code:`INCLUDE_DIRS=...` and/or :code:`LIBRARY_DIRS=...`.

Support
=======
If you are having problems, write to support@geomodelr.com.

License
=======
This project is licensed under the Affero GPL license https://www.gnu.org/licenses/


