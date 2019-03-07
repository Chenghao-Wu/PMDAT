#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

from SMMSAT.src.cython_func.block_loop import *
from SMMSAT.src.cython_func.species_loop import *
from SMMSAT.src.method.analysis import *

class mean_squared_internal_distance(Analysis):

    def __init__(self,List,filename):
        self.List=List
        if not List.AtomListType ==  "species" and not List.AtomListType ==  "type_species":
            print("ERROR::mean_squared_internal_distance: Please choose ListType = type_species")
            return None
        self.filename=filename
        self.SpeciesName=self.List.SpeciesName

        #constructor:
        self.NumberSpeciesAtoms=np.sum(self.List.System.sys_SpeciesDict[self.SpeciesName].AtomsList)
        self.NumberSpecies=self.List.System.sys_SpeciesDict[self.SpeciesName].NumberSpecies

        self.MSID=np.zeros(self.NumberSpeciesAtoms,dtype=np.float32)
        self.weighting=np.zeros(self.List.System.get_NumberSteps)
        print("\nInitializing Mean Squared Internal Distance "+str(self.SpeciesName))

    def excute_Analysis(self):
        print("\nCalculating mean squared internal distance.\n")
        start=time.time()
        block_loop(self,self.List.unwrap_pos)
        self.postprocess_list()
        self.write()
        print("Writing mean squared internal distance to file "+self.filename)
        end=time.time()
        print("\nCalculated mean squared internal distance in " +"{0:.2f}".format(end-start) +" seconds.")
        

    def list_statickernel(self,time,current_trajectory):
        totaldistance=np.sum(SpeciesLoop(self,current_trajectory),axis=0)
        for beadii in range(len(self.MSID)):
            self.MSID[beadii]+=np.mean(totaldistance.diagonal(offset=beadii))
        self.weighting[time]+=1
        
    def listkernel(self,atomii,species_trajetcory):
        distance_atom=np.linalg.norm(species_trajetcory[atomii]-species_trajetcory,axis=1)
        return distance_atom

    def postprocess_list(self):
        self.MSID=self.MSID/np.sum(self.weighting)/self.NumberSpecies
        for index in range(1,len(self.MSID)):
            self.MSID[index]=np.power(self.MSID[index],2)/index
        

    def write(self):
        end_end_distance_file=open(self.filename+".csv","w")
        with end_end_distance_file:
            writer=csv.writer(end_end_distance_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["mean squared internal distance data created by SMMSAT "+self.List.System.Version])
            for index,MSID in enumerate(self.MSID):
                writer.writerow(["{0:.4f}".format(index),("{0:."+str(math_tool.data_precision)+"f}").format(self.MSID[index])])

    