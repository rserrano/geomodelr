"""
	Geomodelr query tool. Tools for using a geomodelr.com geological model.
	Copyright (C) 2016 Geomodelr, Inc.
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.
	
	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
from setuptools import setup, Extension
import os

def if_env(**kwargs):
    varname, default = kwargs.popitem()
    try:
        values = os.environ[varname].split(os.pathsep)
        if len(values) < 1:
           return  default
        return values
    except KeyError:
        return default

def_libraries    = ['boost_python']
def_library_dirs = ['/usr/local/lib', '/usr/lib/x86_64-linux-gnu/']
def_include_dirs = ['/usr/local/include']

cppextension = Extension("geomodelr.cpp", [ 'geomodelr/cpp/match.cpp', 'geomodelr/cpp/geomodel.cpp', 'geomodelr/cpp/model.cpp', 'geomodelr/cpp/section.cpp',], 
                         depends=['geomodelr/cpp/geomodel.hpp'],
                         include_dirs=if_env(INCLUDE_DIRS=def_include_dirs),
                         library_dirs=if_env(LIBRARY_DIRS=def_library_dirs), 
                         libraries=if_env(LIBRARIES=def_libraries),
                         extra_compile_args=['-std=c++11'])

setup(name='geomodelr',
    version='0.1.3',
    description='Geomodelr is the open source query tool for geomodelr.com models.',
    url='http://github.com/rserrano/geomodelr',
    author='Ricardo Serrano',
    author_email='rserrano@geomodelr.com',
    license='AGPL',
    packages=['geomodelr'],
    ext_modules=[cppextension],
    install_requires=['numpy', 'scipy', 'shapely'],
    keywords=['geology', 'geological modelling', 'cross sections', 'geomodelr'],
    entry_points = {
        'console_scripts': [
            'geomodelr=geomodelr.__main__:main'
        ]
    },
    zip_safe=False)

