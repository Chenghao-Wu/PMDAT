#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

from SMMSAT.cython_func.displacement_loop import *
from SMMSAT.cython_func.distance import *
from analysis.analysis import *

class mean_squared_displacement(Analysis):

    def __init__(self,List,filename):
        self.List=List
        self.filename=filename
        self.msd=np.zeros(self.List.System.c_NumberTimeGaps,dtype=np.float64)
        self.weighting=np.zeros(self.List.System.c_NumberTimeGaps)
        print("\nMean Squared displacement")

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
        squeared_distance=calculte_squared_distance(thisii_array,nextii_array)
        self.weighting[timegap]+=1
        self.msd[timegap]+=squeared_distance

    def postprocess_list(self):
        for timeii in range(self.List.System.c_NumberTimeGaps):
            self.msd[timeii]/=float(self.weighting[timeii])

    def write(self):
        MSD_file=open(self.filename+".csv","w")
        with MSD_file:
            writer=csv.writer(MSD_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["MSD data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.c_TimeGap):
                writer.writerow(["{0:.4f}".format(time),("{0:."+str(math_tool.data_precision)+"f}").format(self.msd[index])])


    