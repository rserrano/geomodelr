from setuptools import setup, Extension

cppextension = Extension("cpp", ['geomodelr/cpp/model.cpp'], 
                         include_dirs=['/usr/local/include'], 
                         library_dirs=[ '/usr/local/lib', '/usr/lib/x86_64-linux-gnu/' ], 
                         libraries=['boost_python'])

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

