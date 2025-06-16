from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
from glob import glob
import re, sys

with open('orient/_version.py', 'r') as f:
    __version__ = re.match(r"""__version__ = ["'](.*?)['"]""", f.read()).group(1)


libraries = ['gsl', 'gslcblas']
if sys.platform != 'darwin': libraries.append('stdc++fs')

ext_modules = [
    Pybind11Extension(
        'orient.orient',
        glob('orient/src/*.cpp'),
        include_dirs=['orient/src'],
        cxx_std=14,
        libraries=libraries
    ),
]

setup(
    name='orient',
    version=__version__,
    ext_modules=ext_modules,
    packages=['orient'],
    package_data={
        '': ['data/*.dat']
    }
)
