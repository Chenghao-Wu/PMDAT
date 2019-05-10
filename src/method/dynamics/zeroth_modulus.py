#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
import csv

from SMMSAT.src.cython_func import math_tool

"""
Citation: Paul Müller (2012) Python multiple-tau algorithm
"""
import SMMSAT.multipletau as mpt
from SMMSAT.src.method.analysis import *
class zeroth_modulus(Analysis):
    def __init__(self,List,filename):
        self.List=List
        self.filename=filename
        
    def excute_Analysis(self):

        print("\nCalculating zeroth modulus.\n")
        start=time.time()
        self.kernel()
        self.write()
        print("Writing msd to file "+self.filename)
        end=time.time()
        print("\nCalculated zeroth modulus in " +"{0:.2f}".format(end-start) +" seconds.")

    def kernel(self):
        Txy=self.List.System.sys_data["x"].groupby(level=[0,1]).sum()
        Txz=self.List.System.sys_data["y"].groupby(level=[0,1]).sum()
        Tyz=self.List.System.sys_data["z"].groupby(level=[0,1]).sum()
        # Set m-value for correlation.
        # Gives the number of entries for each inverval of x of 2^x.
        mval = 16

        # Perform autocorrelation. m gives the number of entries for each inverval of x of 2^x.
        self.ACF_xy = mpt.autocorrelate(Txy, m=mval, normalize=False)
        self.ACF_xz = mpt.autocorrelate(Txz, m=mval, normalize=False)
        self.ACF_yz = mpt.autocorrelate(Tyz, m=mval, normalize=False)
    
        self.ACF_av = (self.ACF_xy[:,1]+self.ACF_xz[:,1]+self.ACF_yz[:,1])/3

        one = mpt.autocorrelate(np.ones(Txy.shape[0]), m=mval, normalize=False)

        self.ACF_xy[:,1] = self.ACF_xy[:,1]/one[:,1]
        self.ACF_xz[:,1] = self.ACF_xz[:,1]/one[:,1]
        self.ACF_yz[:,1] = self.ACF_yz[:,1]/one[:,1]
        
        self.ACF_av[:] = self.ACF_av[:]/one[:,1]

    def write(self):
        ZM_file=open(self.filename+".csv","w")
        with ZM_file:
            writer=csv.writer(ZM_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Zeroth Modulus data created by SMMSAT "+self.List.System.Version])
            for index in range(self.ACF_xy.shape[0]):
                writer.writerow(["{0:.4f}".format(self.ACF_xy[index,0]*self.List.System.sys_TimeUnit),("{0:."+str(math_tool.data_precision)+"f}").format(self.ACF_av[index])])

    