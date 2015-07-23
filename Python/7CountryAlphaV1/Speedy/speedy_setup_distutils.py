from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy import get_include
from os import system

shared_obj="speedy.c -fPIC -c -o speedy.o"
print shared_obj
system(shared_obj)

ext_modules=[Extension(
                                    "cython_speedy",
                                    ["cython_speedy.pyx"],

                                    extra_link_args=["speedy.o"])]

setup(name="cython_speedy",
        cmdclass={'build_ext':build_ext},
        include_dirs=[get_include()],
        ext_modules=ext_modules)
