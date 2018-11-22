#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
import csv

from math_func import math_tool

class Correlation_2d(object):
    def __init__(self):
        self.timegap_weighting=np.zeros(0)
        self.firsttime=0
        self.lasttime=0
        self.n_atoms_represented=0
        self.correlation=np.zeros(0)
        self.filename=str(0)
        self.timetable=np.zeros(0)
        self.weighting=np.zeros(self.List.System.c_NumberTimeGaps)
    def excute(self):

        self.excute_Analysis()

    def list_displacementkernel(self):
        pass
    
    def postprocess_list(self):
        for timeii in range(self.firsttime,self.lasttime+1):
            #   print(self.correlation[timeii])
            self.correlation[timeii]/=(float(self.n_atoms_represented)*self.weighting[timeii])
            #print(self.correlation[timeii])

    