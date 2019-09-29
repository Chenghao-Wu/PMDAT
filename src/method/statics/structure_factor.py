#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import time

from SMMSAT.src.cython_func.block_loop import *
from SMMSAT.src.cython_func import math_tool
from SMMSAT.src.method.analysis import *
from SMMSAT.src.method.wave_density import * 


class structure_factor(Analysis):
    def __init__(self,List,filename,MaxLengthScale,geometry,inner=0,outer=-1,fulltrj=False):
        self.List=List
        self.wave_density=WaveDensity(self.List,MaxLengthScale,geometry,inner=0,outer=-1)
        self.structure_factor=np.zeros(self.wave_density.n_wavenumbers,dtype=np.float64)
        block_loop(self.wave_density,self.List.unwrap_pos)
        self.fulltrj=fulltrj
        #print(self.wave_density.density[0])
        self.filename=filename

    def excute_Analysis(self):
        print("\nCalculating radial distribution function.\n")
        start=time.time()
        
        for currenttime in range(self.List.System.get_NumberBlocks):
            self.list_statickernel(currenttime)


        self.postprocess_list()
        self.write()
        print("Writing rdf to file "+self.filename)
        end=time.time()
        print("\nCalculated radial distribution function in " +"{0:.2f}".format(end-start) +" seconds.")
    
    def list_statickernel(self,currenttime):
        for wavenumberii in range(self.wave_density.n_wavenumbers):
            vectorcount = self.wave_density.WaveVector.vectorcount(wavenumberii)
            for vectorii in range(vectorcount):
                real=self.wave_density.density[currenttime][wavenumberii][vectorii].real
                image=self.wave_density.density[currenttime][wavenumberii][vectorii].imag

                self.structure_factor[wavenumberii]+=(real*real+image*image)/vectorcount
    
    def postprocess_list(self):
        self.structure_factor=self.structure_factor/self.List.System.get_NumberAtoms/self.List.System.get_NumberBlocks

    def write(self):
        structurefactor_file=open(self.filename+".csv","w")
        with structurefactor_file:
            writer=csv.writer(structurefactor_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["structure factor by SMMSAT "+self.List.System.Version])
            for wavenumberii in range(self.wave_density.n_wavenumbers):
                writer.writerow(["{0:.4f}".format(self.wave_density.WaveVector.show_ApproxWaveNumber(wavenumberii)),("{0:."+str(math_tool.data_precision)+"f}").format(self.structure_factor[wavenumberii])])
