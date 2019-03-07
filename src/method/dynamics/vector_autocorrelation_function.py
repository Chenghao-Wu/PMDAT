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

import sys
#need revised
class vector_autocorrelation_function(Correlation_2d):

    def __init__(self,List,filename,plane):
        self.List=List
        #check List type
        if List.ListType=="MultibodyList":
            print("\nReading Multibody List...")
        else:
            print("\nERROR:vector_autocorrelation_function:: Reading Wrong Type of List, Please Use Multibody List")
            sys.exit(2)

        self.filename=filename
        self.plane=plane

        self.firsttime = 0
        self.lasttime = self.List.System.get_NumberTimeGaps-1
        self.timetable=self.List.System.get_TimeGapTable

        #construction
        self.weighting=np.zeros(self.List.System.get_NumberTimeGaps)
        self.correlation = np.zeros(self.List.System.get_NumberTimeGaps,dtype=np.float32)
        
        self.timetable=self.List.System.get_TimeGapTable
        self.n_atoms_represented=self.List.sys_data.shape[0]/self.List.System.sys_Reader.sys_NumberFrames
        print("\nInitializing Vector Autocorrelation Function "+plane)

    def excute_Analysis(self):
        print("\nCalculating Vector Autocorrelation Function.\n")
        start=time.time()
        displacement_loop(self,self.List.sys_data[["vector_x","vector_y","vector_z"]])
        self.postprocess_list()
        self.write()
        print("Writing baf to file "+self.filename)
        end=time.time()
        print("\nCalculated Vector Autocorrelation Function in " +"{0:.2f}".format(end-start) +" seconds.")
        
    def list_displacementkernel(self,timegap,thisii_array,nextii_array):
        tempcorrelation=0.0
        if self.plane == "xyz":
            tempcorrelation=np.sum(np.sum(np.multiply(thisii_array,nextii_array),axis=1)/(np.linalg.norm(thisii_array,axis=1)*np.linalg.norm(nextii_array,axis=1)))
            temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(thisii_array),math_tool.unit_vector(nextii_array)),axis=1),2))
        elif self.plane == "xy":
            tempcorrelation=np.sum(np.sum(np.multiply(thisii_array[:,[0,1]],nextii_array[:,[0,1]]),axis=1)/(np.linalg.norm(thisii_array[:,[0,1]],axis=1)*np.linalg.norm(nextii_array[:,[0,1]],axis=1)))
            temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(thisii_array[:,[0,1]]),math_tool.unit_vector(nextii_array[:,[0,1]])),axis=1),2))
        elif self.plane == "xz":
            tempcorrelation=np.sum(np.sum(np.multiply(thisii_array[:,[0,2]],nextii_array[:,[0,2]]),axis=1)/(np.linalg.norm(thisii_array[:,[0,2]],axis=1)*np.linalg.norm(nextii_array[:,[0,2]],axis=1)))
            temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(thisii_array[:,[0,2]]),math_tool.unit_vector(nextii_array[:,[0,2]])),axis=1),2))
        elif self.plane == "yz":
            tempcorrelation=np.sum(np.sum(np.multiply(thisii_array[:,[1,2]],nextii_array[:,[1,2]]),axis=1)/(np.linalg.norm(thisii_array[:,[1,2]],axis=1)*np.linalg.norm(nextii_array[:,[1,2]],axis=1)))
            temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(thisii_array[:,[1,2]]),math_tool.unit_vector(nextii_array[:,[1,2]])),axis=1),2))
        else:
            print("\nERROR:vector_autocorrelation_function::list_displacementkernel, please set right plane")
        self.correlation[timegap] += float(temporientationalcorrelation)
        self.weighting[timegap]+=1

    def write(self):
        CORRELATION_file=open(self.filename+".csv","w")
        with CORRELATION_file:
            writer=csv.writer(CORRELATION_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["correlation data created by SMMSAT "+self.List.System.Version])
            for timeii in range(self.firsttime,self.lasttime+1):
                writer.writerow(["{0:.4f}".format(self.timetable[timeii]),("{0:."+str(math_tool.data_precision)+"f}").format(self.correlation[timeii])])
    