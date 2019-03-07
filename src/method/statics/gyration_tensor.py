#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""


from SMMSAT.src.method.analysis import *
from SMMSAT.src.cython_func.block_loop import *
from SMMSAT.src.cython_func.species_loop import *
np.set_printoptions(precision=math_tool.data_precision)


class gyration_tensor(Analysis):

    def __init__(self,List,filename):
        self.List=List
        if not List.AtomListType ==  "species" and not List.AtomListType ==  "type_species":
            print("ERROR::mean_squared_internal_distance: Please choose ListType = type_species")
            return None

        self.filename=filename
        self.SpeciesName=self.List.SpeciesName

        self.NumberSpeciesAtoms=np.sum(self.List.System.sys_SpeciesDict[self.SpeciesName].AtomsList)
        self.NumberSpecies=self.List.System.sys_SpeciesDict[self.SpeciesName].NumberSpecies

        self.weighting=np.zeros(self.List.System.get_NumberBlocks)
        self.gyrationtensor=np.zeros(self.List.System.get_NumberBlocks,dtype=np.ndarray)
        print("\nInitializing Gyration Tensor "+str(self.SpeciesName))

    def excute_Analysis(self):
        print("\nCalculating gyration tensor.\n")
        start=time.time()
        block_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing gyration tensor to file "+self.filename)
        end=time.time()
        print("\nCalculated gyration tensor in " +"{0:.2f}".format(end-start) +" seconds.")
        

    def list_statickernel(self,time,current_trajectory):
        result=np.sum(np.sum(SpeciesLoop(self,current_trajectory),axis=0),axis=0)
        self.weighting[time]+=self.List.System.sys_SpeciesDict[self.SpeciesName].NumberSpecies
        
        self.gyrationtensor[time]+=result/(np.power(self.NumberSpeciesAtoms,2)*2)

    def listkernel(self,atomii,Trj):
        S=np.dot((Trj[atomii]-Trj).transpose(),(Trj[atomii]-Trj))
        return S

    def postprocess_list(self):
        for timeii in range(self.List.System.get_NumberBlocks):
            self.gyrationtensor[timeii]/=(float(self.weighting[timeii]))
        
    def write(self):
        gyrationtensor_file=open(self.filename+".csv","w")
        with gyrationtensor_file:
            writer=csv.writer(gyrationtensor_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["gyration tensor data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.get_TimeGapTable):
                """
                xx, yy, zz, xy, xz, yz
                """
                writer.writerow(["{0:.4f}".format(time),self.gyrationtensor[index][0][0],self.gyrationtensor[index][1][1],self.gyrationtensor[index][2][2],self.gyrationtensor[index][0][1],self.gyrationtensor[index][0][2],self.gyrationtensor[index][1][2]])


    