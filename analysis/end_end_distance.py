#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

from SMMSAT.cython_func.block_loop import *
from SMMSAT.cython_func.listloop import *
from analysis.analysis import *

class end_end_distance(Analysis):

    def __init__(self,List,filename):
        self.List=List
        self.filename=filename
        self.weighting=np.zeros(self.List.System.c_NumberSteps)
        self.EndEndDistance=np.zeros(self.List.System.c_NumberSteps,dtype=np.float64)

    def excute_Analysis(self):
        print("\nCalculating end to end distance.\n")
        start=time.time()
        timeseris_loop(self,self.List.vector)
        self.postprocess_list()
        self.write()
        print("Writing end to end distance to file "+self.filename)
        end=time.time()
        print("\nCalculating end end distance in " +"{0:.2f}".format(end-start) +" seconds.")
        

    def list_statickernel(self,time,current_trajectory):
        distance=np.linalg.norm(current_trajectory,axis=1)
        self.weighting[time]+=1
        self.EndEndDistance[time]+=np.mean(distance)

    def postprocess_list(self):
        for timeii in range(self.List.System.c_NumberTimeGaps):
            self.EndEndDistance[timeii]/=float(self.weighting[timeii])

    def write(self):
        end_end_distance_file=open(self.filename+".csv","w")
        with end_end_distance_file:
            writer=csv.writer(end_end_distance_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["end_end_distance data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.c_TimeList):
                writer.writerow(["{0:.4f}".format(time),("{0:."+str(math_tool.data_precision)+"f}").format(self.EndEndDistance[index])])


    