#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import time

from SMMSAT.src.cython_func.displacement_loop import *
from SMMSAT.src.cython_func.distance import *
from SMMSAT.src.method.wave_vectors import *
from SMMSAT.src.method.dynamics.correlation_2d import *
from SMMSAT.src.cython_func import math_tool

class intermediate_scattering_function(Correlation_2d):

    def __init__(self,List,filename,index,geometry,MaxLengthScale):
        self.List=List
        self.filename=filename
        self.index=index
        self.geometry=geometry
        self.weighting=np.zeros(self.List.System.get_NumberTimeGaps)
        self.WaveVector=WaveVector(self.List.System,geometry,MaxLengthScale)
        self.correlation = np.zeros(self.List.System.get_NumberTimeGaps,dtype=np.float32)
        self.firsttime = 0
        self.lasttime = self.List.System.get_NumberTimeGaps-1
        self.timetable=self.List.System.get_TimeGapTable
        self.n_atoms_represented=self.List.unwrap_pos.shape[0]/self.List.System.get_NumberFrames
        print("\nInitializing Intermediate Scattering Function "+str(index)+" "+geometry+" "+str(MaxLengthScale))

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
        if self.geometry=="xyz":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array-thisii_array,self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="xy":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[0,1]]-thisii_array[:,[0,1]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="yz":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[1,2]]-thisii_array[:,[1,2]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="xz":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[0,2]]-thisii_array[:,[0,2]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="x":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[0]]-thisii_array[:,[0]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="y":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[1]]-thisii_array[:,[1]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
        elif self.geometry=="z":
            for vectorii in range(self.WaveVector.vectorcount(self.index)):
                tempcorrelation+=np.sum(np.cos(np.dot(nextii_array[:,[2]]-thisii_array[:,[2]],self.WaveVector.vectorlist(self.index).loc[vectorii].values)))/self.WaveVector.vectorcount(self.index)
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

    