#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:24:59 2018

@author: bruce
"""

import cython
cimport cython

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free


def block_loop(analysis,Trj):
    cdef int blockii
    cdef int thisii
    cdef int n_exponentials=analysis.List.System.get_NumberBlocks
    cdef int n_exponential_steps=analysis.List.System.get_BlockSize
    cdef np.ndarray thisii_array
    for blockii in range(n_exponentials):
        thisii = n_exponential_steps*blockii
        thisii_array=Trj.loc[thisii].values
        analysis.list_statickernel(blockii,thisii_array)
        
  