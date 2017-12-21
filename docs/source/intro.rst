
Geomodelr, Introduction.
************************
 
To use geomodelr query tool you just need to::

    import geomodelr
    # load your model.
    model = geomodelr.model_from_file('/path/to/your/model_version.json')
    # query your model.
    unit, gmlr_distance = model.closest((1000, 1000, 0.0))
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

The requirements of Geomodelr Query Tool are:
- Currently, Linux and Mac OS X are supported in Python 2.7 but we plan to support Windows and Python 3.5 in the near future.
- C++ Build tools that support C++11.
- Boost Libraries.
- numpy, scipy and (pip will install them).

In general, you can install geomodelr by calling::

    pip install geomodelr

Ubuntu Linux
------------
You can install it from the command line::

    # This will install boost.
    sudo apt-get install libboost-all-dev
    # This will install geomodelr globally. 
    sudo pip install geomodelr
    # OR You can also use virtualenv.
    virtualenv env && source env/bin/activate
    pip install geomodelr

Mac OS X
--------
- Install mac ports from https://www.macports.org/
- Now from the command line::

    sudo ports install boost
    sudo ports install pip
    INCLUDE_DIRS="/opt/local/include/boost" LIBRARY_DIRS="/opt/local/lib" LIBRARIES="boost_python-mt" pip install geomodelr --user

Mac OS X El Capitan has a binary wheel so you don't need to install boost or anything besides pip and geomodelr. 

Support
=======
If you are having problems, write to support@geomodelr.com.

License
=======
This project is licensed under the Affero GPL license https://www.gnu.org/licenses/


