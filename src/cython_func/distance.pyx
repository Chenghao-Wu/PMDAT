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

def calculate_distance(np.ndarray trajectory_1, np.ndarray trajectory_2):
    cdef int rows=trajectory_1.shape[0]
    cdef int columns=trajectory_1.shape[1]
    cdef np.ndarray distance_array=np.zeros([rows,1],dytype=np.float64)
    distance_array=np.sqrt(np.sum(np.power(trajectory_1-trajectory_2,2),axis=1),axis=1)
    cdef float distance
    distance=np.sum(distance_array,axis=0)[0]/rows
    return distance

def calculte_squared_distance(np.ndarray trajectory_1, np.ndarray trajectory_2):
    cdef int rows=trajectory_1.shape[0]
    cdef int columns=trajectory_1.shape[1]
    cdef np.ndarray distance_array=np.zeros([rows,1],dtype=np.float64)
    distance_array=np.sum(np.power(trajectory_1-trajectory_2,2),axis=1)
    cdef float distance
    distance=np.sum(distance_array,axis=0)/rows
    return distance
