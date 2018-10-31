#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import time

from SMMSAT.cython_func.block_loop import *
from analysis.analysis import *

#need revised
class order_parameter(Analysis):

    def __init__(self,List,filename,plane,n_bins,fullblock=0):
        self.List=List
        self.filename=filename
        self.plane=plane
        self.n_bins=n_bins
        self.fullblock=fullblock
        self.weighting=np.zeros(self.List.System.c_NumberBlocks)
        self.order_parameter = np.zeros(self.List.System.c_NumberBlocks,dtype=np.ndarray)
        for blockii in range(self.List.System.c_NumberBlocks):
            self.order_parameter[blockii]=np.zeros(self.n_bins)
        if self.plane == "x":
            self.MaxDistance=self.List.System.c_AbsBoxSize[0]
            self.corr_vector=np.array([1,0,0])
        elif self.plane == "y":
            self.MaxDistance=self.List.System.c_AbsBoxSize[1]
            self.corr_vector=np.array([0,1,0])
        elif self.plane == "z":
            self.MaxDistance=self.List.System.c_AbsBoxSize[1]
            self.corr_vector=np.array([0,0,1])
        else:
            print("\nERROR:Order_Parameter::__init__, please set plane in 'x', 'y' or 'z'")
        
        self.dr=self.MaxDistance/self.n_bins
        self.hist_bins=[-self.MaxDistance/2+self.dr*i for i in range(n_bins+1)]
        print("\nOrder Parameter "+plane+" "+str(n_bins))

    def excute_Analysis(self):
        print("\nCalculating Order_Parameter.\n")
        start=time.time()
        block_loop(self,self.List.vector)
        self.postprocess_list()
        self.write()
        print("Writing msd to file "+self.filename)
        end=time.time()
        print("\nCalculated Order_Parameter in " +"{0:.2f}".format(end-start) +" seconds.")
        
    def list_statickernel(self,block,current_trajectry):
        self.weighting[block]+=1
        for index_bin in range(self.n_bins):
            start_bin=self.hist_bins[index_bin]
            end_bin=self.hist_bins[index_bin+1]
            if self.plane == "x":
                bin_trj=current_trajectry[(self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]>start_bin) & (self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]<end_bin)]
                tempcorrelation=np.sum(np.multiply(bin_trj,self.corr_vector),axis=1)/(np.linalg.norm(bin_trj,axis=1))
                temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(bin_trj),self.corr_vector),axis=1),2))/bin_trj.shape[0]
            elif self.plane == "y":
                bin_trj=current_trajectry[(self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]>start_bin) & (self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]<end_bin)]
                tempcorrelation=np.sum(np.multiply(bin_trj,self.corr_vector),axis=1)/(np.linalg.norm(bin_trj,axis=1))
                temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(bin_trj),self.corr_vector),axis=1),2))/bin_trj.shape[0]
            elif self.plane == "z":
                bin_trj=current_trajectry[(self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]>start_bin) & (self.List.wrap_pos.loc[block*self.List.System.c_BlockSize].values[:,0]<end_bin)]
                tempcorrelation=np.sum(np.multiply(bin_trj,self.corr_vector),axis=1)/(np.linalg.norm(bin_trj,axis=1))
                temporientationalcorrelation=np.sum(math_tool.legendre_polynomial(np.sum(np.multiply(math_tool.unit_vector(bin_trj),self.corr_vector),axis=1),2))/bin_trj.shape[0]
            else:
                print("\nERROR:Order_Parameter::list_displacementkernel, please set right plane")
            self.order_parameter[block][index_bin] += temporientationalcorrelation
        

    def postprocess_list(self):
        self.order_parameter=np.sum(self.order_parameter,axis=0)/np.sum(self.weighting)

    def write(self):
        CORRELATION_file=open(self.filename+".csv","w")
        with CORRELATION_file:
            writer=csv.writer(CORRELATION_file,delimiter='\t',
                            quotechar='\t', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["correlation data created by SMMSAT "+self.List.System.Version])
            for binii in range(self.n_bins):
                writer.writerow(["{0:.4f}".format(self.hist_bins[binii]),("{0:."+str(math_tool.data_precision)+"f}").format(self.order_parameter[binii])])
    