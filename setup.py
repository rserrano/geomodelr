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

cppextension = Extension("cpp", ['geomodelr/cpp/model.cpp', 'geomodelr/cpp/section.cpp'], 
                         include_dirs=['/usr/local/include'], 
                         library_dirs=[ '/usr/local/lib', '/usr/lib/x86_64-linux-gnu/' ], 
                         libraries=['boost_python'], 
                         extra_compile_args=['-std=c++11'])

setup(name='geomodelr',
      version='0.1',
      description='GeoModelR is the open source client of the geological modeler.',
      url='http://github.com/rserrano/geomodelr',
      author='Ricardo Serrano',
      author_email='rserrano@geomodelr.com',
      license='AGPL',
      packages=['geomodelr'],
      ext_modules=[cppextension],
      install_requires=['sortedcontainers', 'numpy', 'scipy', 'shapely'],
      zip_safe=False)

