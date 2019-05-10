#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import numpy as np

from SMMSAT.src.method.analysis import *
from SMMSAT.src.cython_func.block_loop import *
import sys
class end_end_distance(Analysis):

    def __init__(self,List,filename,calc_Dist=False):
        self.List=List
        if not List.AtomListType ==  "species" and not List.AtomListType ==  "type_species":
            print("ERROR::end_end_distance: Please choose ListType = type_species")
            sys.exit(2)
        self.filename=filename
        self.SpeciesName=self.List.SpeciesName
        self.NumberSpecies=self.List.System.sys_SpeciesDict[self.SpeciesName].NumberSpecies
        self.SpeciesLength=self.List.selectedspecieslength
        self.SpeciesIndex=list(self.List.System.sys_SpeciesDict.keys()).index(self.SpeciesName)
        
        #constructor:
        self.weighting=np.zeros(self.List.System.get_NumberSteps)
        self.EndEndDistance=np.zeros(self.List.System.get_NumberSteps,dtype=np.float32)
        self.Dist_EndtoEndDistance=np.zeros(self.List.System.get_NumberSteps,dtype=np.ndarray)
        self.Bins_EndtoEndDistance=np.arange(np.sum(self.List.System.sys_SpeciesDict[self.SpeciesName].AtomsList)*1.3+1000*(1/(np.sum(self.List.System.sys_SpeciesDict[self.SpeciesName].AtomsList))))
        
        self.calc_Dist=calc_Dist

        if self.calc_Dist == True: 
            print("\nInitializing End to End Distance with Distribution "+str(self.SpeciesName))
        elif self.calc_Dist == False:
            print("\nInitializing End to End Distance "+str(self.SpeciesName))

    def excute_Analysis(self):
        print("\nCalculating end to end distance.\n")
        start=time.time()
        block_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing end to end distance to file "+self.filename)
        end=time.time()
        print("\nCalculated end end distance in " +"{0:.2f}".format(end-start) +" seconds.")

    def list_statickernel(self, blockii, current_trajectory):
        distance=np.zeros(self.NumberSpecies)
        for monomerii in range(self.NumberSpecies):
            distance_monomer=np.linalg.norm(current_trajectory[self.SpeciesLength*monomerii:self.SpeciesLength+self.SpeciesLength*monomerii][-1]-current_trajectory[self.SpeciesLength*monomerii:self.SpeciesLength+self.SpeciesLength*monomerii][0])
            distance[monomerii]=distance_monomer
        self.EndEndDistance[blockii]+=np.mean(distance)
        if self.calc_Dist==True:
            Dist,Bins=np.histogram(distance,bins=self.Bins_EndtoEndDistance,density=True)
            self.Dist_EndtoEndDistance[blockii]+=Dist
        
    def postprocess_list(self):
        if self.calc_Dist==True:
            self.Dist_EndtoEndDistance=np.sum(self.Dist_EndtoEndDistance,axis=0)
            self.Dist_EndtoEndDistance/=self.List.System.get_NumberBlocks

    def write(self):
        end_end_distance_file=open(self.filename+".csv","w")
        with end_end_distance_file:
            writer=csv.writer(end_end_distance_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["end_end_distance data created by SMMSAT "+self.List.System.Version])
            for index,time in enumerate(self.List.System.get_TimeGapTable):
                writer.writerow(["{0:.4f}".format(time),("{0:."+str(math_tool.data_precision)+"f}").format(self.EndEndDistance[index])])
        if self.calc_Dist==True:
            dist_end_end_distance_file=open(self.filename+"_dist.csv","w")
            with dist_end_end_distance_file:
                writer=csv.writer(dist_end_end_distance_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["end_end_distance distribution data created by SMMSAT "+self.List.System.Version])
                for index in range(len(self.Bins_EndtoEndDistance)-1):
                    writer.writerow([self.Bins_EndtoEndDistance[index],np.round(self.Dist_EndtoEndDistance[index],6)])

    