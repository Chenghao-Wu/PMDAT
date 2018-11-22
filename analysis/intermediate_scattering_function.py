#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import time

from SMMSAT.cython_func.displacement_loop import *
from SMMSAT.cython_func.distance import *
from analysis.wave_vectors import *
from analysis.correlation_2d import *

from math_func import math_tool

class intermediate_scattering_function(Correlation_2d):

    def __init__(self,List,filename,index,plane,MaxLengthScale,fullblock=0):
        self.List=List
        self.filename=filename
        self.index=index
        self.fullblock=fullblock
        self.weighting=np.zeros(self.List.System.c_NumberTimeGaps)
        self.WaveVector=WaveVector(self.List.System,plane,MaxLengthScale)
        self.correlation = np.zeros(self.List.System.c_NumberTimeGaps,dtype=np.float64)
        self.firsttime = 0
        self.lasttime = self.List.System.c_NumberTimeGaps-1
        self.timetable=self.List.System.c_TimeGap
        self.n_atoms_represented=self.List.unwrap_pos.shape[0]/self.List.System.c_DataFrame.c_NumberFrames
        print("\nIntermediate Scattering Function "+str(index)+" "+plane+" "+str(MaxLengthScale)+" "+str(fullblock))

    def excute_Analysis(self):
        print("\nCalculating intermediate scattering function.\n")
        start=time.time()
        displacement_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing msd to file "+self.filename)
        end=time.time()
        print("\nCalculated intermediate scattering function in " +"{0:.2f}".format(end-start) +" seconds.")
        

    def list_displacementkernel(self,timegap,thisii_array,nextii_array):
        tempcorrelation=0
        for vectorii in range(self.WaveVector.vectorcount(self.index)):
            tempcorrelation+=np.sum(np.cos(np.dot(nextii_array-thisii_array,self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        self.correlation[timegap] += float(tempcorrelation)
        self.weighting[timegap]+=1
    
    def write(self):
        CORRELATION_file=open(self.filename+".csv","w")
        with CORRELATION_file:
            writer=csv.writer(CORRELATION_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["correlation data created by SMMSAT "+self.List.System.Version])
            writer.writerow(["{0:.4f}".format(self.WaveVector.show_ApproxWaveNumber(self.index))])
            for timeii in range(self.firsttime,self.lasttime+1):
                writer.writerow(["{0:.4f}".format(self.timetable[timeii]),("{0:."+str(math_tool.data_precision)+"f}").format(self.correlation[timeii])])

    