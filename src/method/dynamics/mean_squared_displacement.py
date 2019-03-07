#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

from SMMSAT.src.cython_func.displacement_loop import *
from SMMSAT.src.cython_func.distance import *
from SMMSAT.src.cython_func import math_tool
from SMMSAT.src.method.analysis import *


class mean_squared_displacement(Analysis):

    def __init__(self,List,filename,geometry):
        self.List=List
        self.filename=filename
        self.geometry=geometry
        self.MSD=np.zeros(self.List.System.get_NumberTimeGaps,dtype=np.float32)
        self.weighting=np.zeros(self.List.System.get_NumberTimeGaps)
        print("\nInitializing Mean Squared displacement")

    def excute_Analysis(self):
        print("\nCalculating mean square displacement.\n")
        start=time.time()
        displacement_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing msd to file "+self.filename)
        end=time.time()
        print("\nCalculated mean square displacement in " +"{0:.2f}".format(end-start) +" seconds.")
        

    def list_displacementkernel(self,timegap,thisii_array,nextii_array):
        if self.geometry=="xyz":
            squeared_distance=calculte_squared_distance(thisii_array,nextii_array)
        elif self.geometry=="xy":
            squeared_distance=calculte_squared_distance(thisii_array[:,[0,1]],nextii_array[:,[0,1]])
        elif self.geometry=="yz":
            squeared_distance=calculte_squared_distance(thisii_array[:,[1,2]],nextii_array[:,[1,2]])
        elif self.geometry=="xz":
            squeared_distance=calculte_squared_distance(thisii_array[:,[0,2]],nextii_array[:,[0,2]])
        elif self.geometry=="x":
            squeared_distance=calculte_squared_distance(thisii_array[:,[0]],nextii_array[:,[0]])
        elif self.geometry=="y":
            squeared_distance=calculte_squared_distance(thisii_array[:,[1]],nextii_array[:,[1]])
        elif self.geometry=="z":
            squeared_distance=calculte_squared_distance(thisii_array[:,[2]],nextii_array[:,[2]])
        self.weighting[timegap]+=1
        self.MSD[timegap]+=squeared_distance

    def postprocess_list(self):
        for timeii in range(self.List.System.sys_NumberTimeGaps):
            self.MSD[timeii]/=float(self.weighting[timeii])

    def write(self):
        MSD_file=open(self.filename+".csv","w")
        with MSD_file:
            writer=csv.writer(MSD_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["MSD data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.get_TimeGapTable):
                writer.writerow(["{0:.4f}".format(time),("{0:."+str(math_tool.data_precision)+"f}").format(self.MSD[index])])


    