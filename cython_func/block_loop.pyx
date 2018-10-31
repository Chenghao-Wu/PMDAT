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


def block_loop(analysis,Trj):
    cdef int blockii
    cdef int thisii
    cdef int n_exponentials=analysis.List.System.c_NumberBlocks
    cdef int n_exponential_steps=analysis.List.System.c_BlockSize
    cdef np.ndarray thisii_array
    for blockii in range(n_exponentials):
        thisii = n_exponential_steps*blockii
        thisii_array=Trj.loc[thisii].values
        analysis.list_statickernel(blockii,thisii_array)

def timeseris_loop(analysis,Trj):
    cdef int blockii
    cdef int thisii
    cdef int expstepii
    cdef int n_exponentials=analysis.List.System.c_NumberBlocks
    cdef int n_exponential_steps=analysis.List.System.c_BlockSize
    cdef np.ndarray thisii_array
    for blockii in range(n_exponentials):
        for expstepii in range(n_exponential_steps):
            thisii = expstepii+n_exponential_steps*blockii
            thisii_array=Trj.loc[thisii].values
            analysis.list_statickernel(thisii,thisii_array)
    thisii = n_exponentials*n_exponential_steps
    thisii_array=Trj.loc[thisii].values
    analysis.list_statickernel(thisii,thisii_array)