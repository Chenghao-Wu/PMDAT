#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""


from SMMSAT.cython_func.block_loop import *
from SMMSAT.cython_func.listloop import *
from analysis.analysis import *

np.set_printoptions(precision=math_tool.data_precision)


class gyration_tensor(Analysis):

    def __init__(self,List,filename):
        self.List=List
        self.filename=filename
        self.weighting=np.zeros(self.List.System.c_NumberSteps)
        self.gyrationtensor=np.zeros(self.List.System.c_NumberSteps,dtype=np.ndarray)
        for stepii in range(self.List.System.c_NumberSteps):
            self.gyrationtensor[stepii]=np.zeros(6)

    def excute_Analysis(self):
        print("\nCalculating gyration tensor.\n")
        start=time.time()
        timeseris_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing gyration tensor to file "+self.filename)
        end=time.time()
        print("\nCalculating gyration tensor in " +"{0:.2f}".format(end-start) +" seconds.")
        
    def list_statickernel(self,time,current_trajectory):
        result=SpeciesLoop(self,current_trajectory)
        result=np.sum(np.sum(result,axis=1)/(np.power(np.sum(self.List.System.SpeciesList[self.List.SpeciesName].AtomsList),2)*2),axis=0)
        self.weighting[time]+=self.List.System.SpeciesList[self.List.SpeciesName].NumberSpecies
        self.gyrationtensor[time]+=result

    def listkernel(self,atomii,Trj):
        Sxx=np.sum(np.multiply((Trj[atomii][0]-Trj[:,0]),(Trj[atomii][0]-Trj[:,0])),axis=0)
        Syy=np.sum(np.multiply((Trj[atomii][1]-Trj[:,1]),(Trj[atomii][1]-Trj[:,1])),axis=0)
        Szz=np.sum(np.multiply((Trj[atomii][2]-Trj[:,2]),(Trj[atomii][2]-Trj[:,2])),axis=0)
        Sxy=np.sum(np.multiply((Trj[atomii][0]-Trj[:,0]),(Trj[atomii][1]-Trj[:,1])),axis=0)
        Sxz=np.sum(np.multiply((Trj[atomii][0]-Trj[:,0]),(Trj[atomii][2]-Trj[:,2])),axis=0)
        Syz=np.sum(np.multiply((Trj[atomii][1]-Trj[:,1]),(Trj[atomii][2]-Trj[:,2])),axis=0)
        S=np.array([Sxx,Syy,Szz,Sxy,Sxz,Syz])
        return S

    def postprocess_list(self):
        for timeii in range(self.List.System.c_NumberSteps):
            self.gyrationtensor[timeii]/=float(self.weighting[timeii])

    def write(self):
        gyrationtensor_file=open(self.filename+".csv","w")
        with gyrationtensor_file:
            writer=csv.writer(gyrationtensor_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["gyration tensor data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.c_TimeList):
                writer.writerow(["{0:.4f}".format(time),self.gyrationtensor[index][0],self.gyrationtensor[index][1],self.gyrationtensor[index][2],self.gyrationtensor[index][3],self.gyrationtensor[index][4],self.gyrationtensor[index][5]])


    