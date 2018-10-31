#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:28:59 2018

@author: bruce
"""

from pathlib import Path
import numpy as np
import pandas as pd

from analysis.analysis import *

class WaveDensity(Analysis):
    def __init__(self,System,WaveVector,inner,outer):
        self.System=System
        self.WaveVector=WaveVector
        self.first_wavenumber_index = inner
        self.last_wavenumber_index = outer
        self.n_wavenumbers = self.last_wavenumber_index - self.first_wavenumber_index + 1
        self.n_times = self.System.show_n_timesteps()
        self.density=np.zeros(self.List.System.c_NumberBlocks,dtype=np.ndarray)
        self.n_atoms = np.zeros(self.n_times)
        for tii in range(self.n_times):
            self.density[tii]=np.zeros(self.n_wavenumbers,dtype=np.ndarray)
            for wavenumberii in range(self.n_wavenumbers):
                self.density[tii][wavenumberii]=np.zeros(self.WaveVector.vectorcount,dtype=np.complex)

    def list_statickernel(self,currenttime,current_trajectory):
        for wavenumberii in range(self.first_wavenumber_index,self.last_wavenumber_index+1):
            vectorlist = self.WaveVector.vectorlist(wavenumberii)
            vectorcount = self.WaveVector.vectorcount(wavenumberii)
            for vectorii in range(vectorcount):
                k_dot_r=np.dot(current_trajectory,vectorlist[vectorii])
                temcomplex=np.sum(np.cos(k_dot_r))+np.sum(np.sin(k_dot_r))*1j
                self.density[currenttime][wavenumberii-self.first_wavenumber_index][vectorii] +=temcomplex

    