#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:28:59 2018

@author: bruce
"""




import os
import cython
cimport cython

import numpy as np
cimport numpy as np

import pandas as pd

#os.environ["qvectorDir"]=os.environ["HOME"]+"/anaconda3/lib/python3.6/site-packages/amdat"


def compute_k_vectors(float max_length_scale,plane):
    wavegrid_spacing = 2.0*np.pi/max_length_scale
    n_wavenumbers=99
    maxrange = wavegrid_spacing * 100.0
    delta_wavenumber = wavegrid_spacing/2.0	
    #approx_wavenumber=[]
    #for wavenumberii in range(n_wavenumbers):
    #    approx_wavenumber.append(wavegrid_spacing+delta_wavenumber*wavenumberii)
    cdef np.ndarray wavevector_=np.zeros([n_wavenumbers,300],dtype=np.float64)
    if plane == "xyz":
        wavevector_=wavegrid_spacing*np.array(read_vectors())
    return wavevector_

def read_vectors():
    wavevector=[]
    n_wavenumbers=99
    for wavenumberii in range(n_wavenumbers):
        buff=".%03d" % (wavenumberii+2)
        filename="qvector"+str(buff)
        vector=pd.read_csv(os.environ["qvectorDir"]+"/qvectors/qvectors3d/"+filename,header=None,sep="\s+")
        wavevector.append(vector.values)
    return wavevector

def compute_k_values(max_length_scale,plane):
    wavegrid_spacing = 2.0*np.pi/max_length_scale
    n_wavenumbers=99
    maxrange = wavegrid_spacing * 100.0
    delta_wavenumber = wavegrid_spacing/2.0	
    cdef np.ndarray wavevector_=np.zeros([n_wavenumbers,300],dtype=np.float64)
    if plane == "xyz":
        wavevector=wavegrid_spacing*np.array(read_vectors())
    wavenumber=np.zeros(len(wavevector))
    cdef int wavevector_pointer
    cdef int size_wavevectors=len(wavevector)
    for wavevector_pointer in range(size_wavevectors):
        mean_k_values=calculte_mean_k_values(wavevector[wavevector_pointer])
        wavenumber[wavevector_pointer]=mean_k_values
    return wavenumber

    
def calculte_mean_k_values(np.ndarray wavevectors):
    cdef float mean_wavenumber
    #print(wavevectors)
    mean_wavenumber=np.sum(np.sqrt(np.sum(np.power(wavevectors,2),axis=1)))/len(wavevectors)
    return mean_wavenumber

def calculate_stdv_k_values(np.ndarray wavevectors, float mean_wavenumber):
    cdef float stdv_wavenumber
    cdef float deviation
    cdef int wavevector_pointer
    cdef int size_wavevectors=len(wavevectors)
    for wavevector_pointer in range(size_wavevectors):
        deviation+=np.power(np.sqrt(np.sum(np.power(wavevectors[wavevector_pointer],2)))-mean_wavenumber,2)
    stdv_wavenumber=np.sqrt(deviation)/size_wavevectors
    return stdv_wavenumber

