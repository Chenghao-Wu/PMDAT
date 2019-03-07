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
from libcpp.vector cimport vector


def pandas_skiplist(np.ndarray skip_FirstList,int PeriodicSize,int n_period):
    cdef vector[int] vect
    cdef int pointer
    cdef int pointer_2
    cdef int length_FirstList=len(skip_FirstList)
    out_put=[]
    for pointer in range(n_period-1):
        out_put.append(skip_FirstList+(PeriodicSize+length_FirstList)*(pointer+1))
    out_put=np.concatenate((skip_FirstList,out_put),axis=None)
    return out_put
