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


def displacement_loop(analysis,Trj):
    cdef int blockii
    cdef int expii
    cdef int thisii

    cdef int nextii
    cdef int timegapii						
    cdef int block_timegapii
    cdef int displacement_count
    cdef int n_exponential_steps=analysis.List.System.get_BlockSize
    cdef int n_exponentials=analysis.List.System.get_NumberBlocks
    cdef int displacement_limit=0
    cdef np.ndarray thisii_array
    cdef np.ndarray nextii_array

    for timegapii in range(n_exponential_steps):
        displacement_count=0
        for blockii in range(n_exponentials):
            thisii = n_exponential_steps*blockii
            nextii = thisii+timegapii
            thisii_array=Trj.loc[thisii].values
            nextii_array=Trj.loc[nextii].values
            analysis.list_displacementkernel(timegapii,thisii_array,nextii_array)
            displacement_count+=1
            if displacement_count == displacement_limit:
                break

    for timegapii in range(n_exponential_steps,analysis.List.System.get_NumberTimeGaps-1):
        displacement_count=0
        block_timegapii = timegapii - n_exponential_steps + 1
        for blockii in range(n_exponentials-block_timegapii):
            thisii = n_exponential_steps*blockii+expii
            nextii = thisii + n_exponential_steps*block_timegapii
            thisii_array=Trj.loc[thisii].values
            nextii_array=Trj.loc[nextii].values
            
            analysis.list_displacementkernel(timegapii,thisii_array,nextii_array)
            displacement_count+=1
            if displacement_count == displacement_limit:
                break
        if displacement_count == displacement_limit:
            break
    thisii_array=Trj.loc[0].values
    nextii_array=Trj.loc[analysis.List.System.get_NumberSteps-1].values
    analysis.list_displacementkernel(analysis.List.System.get_NumberTimeGaps-1,thisii_array,nextii_array)