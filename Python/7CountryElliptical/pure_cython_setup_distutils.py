from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(name="pure_cython", ext_modules=cythonize('pure_cython.pyx'), include_dirs=[numpy.get_include()])
