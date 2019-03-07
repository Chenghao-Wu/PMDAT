#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import cython
cimport cython

import numpy as np
cimport numpy as np

import pandas as pd

@cython.boundscheck(False)
def AtomLoop(analysis,Trj):
    cdef int atomii
    cdef int n_atoms=analysis.List.wrap_pos.shape[0]/analysis.List.System.get_NumberFrames
    result=[]
    for atomii in range(n_atoms):
        listkernel=analysis.listkernel(atomii,Trj)
        result.append(listkernel)
    return result