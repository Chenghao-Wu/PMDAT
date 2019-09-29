#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:28:59 2018

@author: bruce
"""

from pathlib import Path
import numpy as np
import pandas as pd
import os

from SMMSAT.src.method.wave_vectors import *

PackagePath=os.path.dirname(os.path.realpath(__file__))

class WaveDensity(object):
    def __init__(self,List,MaxLengthScale,geometry,inner=0,outer=-1,fulltrj=False):
        self.List=List
        self.fulltrj=fulltrj
        self.WaveVector=WaveVector(self.List.System,geometry,MaxLengthScale)

        self.first_wavenumber_index = inner
        if outer==-1:
            self.last_wavenumber_index = self.WaveVector.NumberWaveNumbers
        else:
            self.last_wavenumber_index = outer
        self.n_wavenumbers = self.last_wavenumber_index - self.first_wavenumber_index
        self.n_times = self.List.System.get_NumberBlocks
        self.density=np.zeros(self.List.System.get_NumberBlocks,dtype=np.ndarray)
        self.n_atoms = np.zeros(self.n_times)
        for tii in range(self.n_times):
            self.density[tii]=np.zeros(self.n_wavenumbers,dtype=np.ndarray)
            for wavenumberii in range(self.n_wavenumbers):
                self.density[tii][wavenumberii]=np.zeros(self.WaveVector.vectorcount(wavenumberii),dtype=np.complex)

    def list_statickernel(self,currenttime,current_trajectory):
        for wavenumberii in range(self.first_wavenumber_index,self.last_wavenumber_index):
            vectorlist = self.WaveVector.vectorlist(wavenumberii)
            vectorcount = self.WaveVector.vectorcount(wavenumberii)
            for vectorii in range(vectorcount):
                k_dot_r=np.dot(current_trajectory,(vectorlist.loc[vectorii]).values)
                temcomplex=np.sum(np.cos(k_dot_r))+np.sum(np.sin(k_dot_r))*1j
                self.density[currenttime][wavenumberii-self.first_wavenumber_index][vectorii] +=temcomplex
