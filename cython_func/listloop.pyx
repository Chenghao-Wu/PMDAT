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

import pandas as pd
@cython.boundscheck(False)
def ListLoop(analysis,Trj):
    cdef int atomii
    cdef int n_atoms=analysis.List.wrap_pos.shape[0]/analysis.List.System.c_DataFrame.c_NumberFrames
    result=[]
    for atomii in range(n_atoms):
        listkernel=analysis.listkernel(atomii,Trj)
        result.append(listkernel)
    return result
@cython.boundscheck(False)
def SpeciesLoop(analysis,Trj):
    cdef Py_ssize_t speciesii
    cdef Py_ssize_t atomii
    cdef int startii
    cdef int endii
    cdef Py_ssize_t NumberSpecies=analysis.List.System.SpeciesList[analysis.List.SpeciesName].NumberSpecies
    cdef Py_ssize_t n_atoms=np.sum(analysis.List.System.SpeciesList[analysis.List.SpeciesName].AtomsList)
    cdef np.ndarray SpeciesTrj=np.zeros([n_atoms,3],dtype=np.float64)
    result=[]
    for speciesii in range(NumberSpecies):
        startii=speciesii*n_atoms
        endii=(speciesii+1)*(n_atoms)
        species_result=[]
        for atomii in range(n_atoms):
            listkernel=analysis.listkernel(atomii,Trj[startii:endii])
            species_result.append(listkernel)
        result.append(species_result)
    return result