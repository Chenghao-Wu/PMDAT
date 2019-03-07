#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:19:37 2018

@author: zwu
"""

#compile command: python setup.py build_ext --inplace

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy


extensions = [
    
    Extension(
        "create_list", ["create_list.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
            
    Extension(
        "displacement_loop", ["displacement_loop.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
            
    Extension(
        "distance", ["distance.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
    Extension(
        "math_tool", ["math_tool.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
    Extension(
        "block_loop", ["block_loop.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
    Extension(
        "species_loop", ["species_loop.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++"),
    Extension(
        "atom_loop", ["atom_loop.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++")
            ]

setup(
    ext_modules=cythonize(extensions)
)   