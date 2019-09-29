#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:28:59 2018

@author: bruce
"""

from pathlib import Path
import numpy as np
import pandas as pd
import os

PackagePath=os.path.dirname(os.path.realpath(__file__))

class WaveVector(object):
    def __init__(self,System,plane,MaxLengthScale):
        self.System=System
        self.plane=plane
        min_size=float(0)
        if MaxLengthScale<=0:
            print("ERROR:wave_vector::__init__, Maximum lengthscale for wavevector decomposition <= 0; defaulting to box size.\n")
            if plane=="xy":
                if self.System.get_AbsBoxSize[0]<self.System.get_AbsBoxSize[1]:
                    min_size=self.System.get_AbsBoxSize[0]
                else:
                    min_size=self.System.get_AbsBoxSize[1]
            elif plane == "xz":
                if self.System.get_AbsBoxSize[0]<self.System.get_AbsBoxSize[2]:
                    min_size=self.System.get_AbsBoxSize[0]
                else:
                    min_size=self.System.get_AbsBoxSize[2]
            elif plane == "yz":
                if self.System.get_AbsBoxSize[1]<self.System.get_AbsBoxSize[2]:
                    min_size=self.System.get_AbsBoxSize[1]
                else:
                    min_size=self.System.get_AbsBoxSize[2]
            elif plane == "xyz":
                min_size=min(self.System.get_AbsBoxSize[0],self.System.get_AbsBoxSize[1],self.System.get_AbsBoxSize[2])
            elif plane == "x":
                min_size=self.System.get_AbsBoxSize[0]
            elif plane == "y":
                min_size=self.System.get_AbsBoxSize[1]
            elif plane == "z":
                min_size=self.System.get_AbsBoxSize[2]
            else:
                print("ERROR:wave_vector::__init__, plane command not recognized. valid inputs are: x, y, z, xy, xz, yz, xyz")
        else:
            if MaxLengthScale > min(self.System.get_AbsBoxSize[0],self.System.get_AbsBoxSize[1],self.System.get_AbsBoxSize[2]):
                print("ERROR:wave_vector::__init__, Maximum length scale larger than smallest box dimension. Consider choosing a different maximum for wavevector decomposition.")
            min_size = MaxLengthScale
        if min_size == 0:
            print("ERROR:wave_vector::__init__, min_size="+str(min_size))
        self.wavegrid_spacing = 2.0*np.pi/min_size	    #define spacing to be given by wavenumber corresponding to smallest system dimension
        self.maxrange = self.wavegrid_spacing * 100.0		#wavevector grid 100 times this length in each direction from the origin
        self.delta_wavenumber = self.wavegrid_spacing/2.0	#define the thickness of the wavenumber bins to be half a gridspacing
        self.NumberWaveNumbers=299
        self.wavevector=[]
        self.approx_wavenumber=[]
        for wavenumberii in range(self.NumberWaveNumbers):
            self.approx_wavenumber.append(self.wavegrid_spacing+self.delta_wavenumber*wavenumberii)
        if plane=="xyz":
            self.read_vectors_3d()
        if plane == "xy" and plane == "xz" and plane == "yz":
            self.read_vectors_2d(plane)
        if plane == "x" and plane == "y" and plane == "z":
            self.read_vectors_1d(plane)

    def read_vectors_1d(self,plane):
        for wavenumberii in range(self.NumberWaveNumbers):
            buff=".%03d" % (wavenumberii+2)
            filename=PackagePath+"/qvectors/qvectors3d/qvector"+str(buff)
            qvector_file = Path(filename)
            if qvector_file.is_file():
                vector=pd.read_csv(filename,header=None,sep="\s+")
                if plane == "x":
                    x_vector=vector.loc[:,0]
                    y_vector=pd.Series(np.zeros(len(x_vector)))
                    z_vector=pd.Series(np.zeros(len(x_vector)))
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
                elif plane == "y":
                    y_vector=vector.loc[:,0]
                    x_vector=pd.Series(np.zeros(len(x_vector)))
                    z_vector=pd.Series(np.zeros(len(x_vector)))
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
                elif plane == "z":
                    z_vector=vector.loc[:,0]
                    x_vector=pd.Series(np.zeros(len(x_vector)))
                    y_vector=pd.Series(np.zeros(len(x_vector)))
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
            else:
                print("ERROR:wave_vector::read_vectors_1d, "+filename+" can not open")
                self.wavevector.append(vector)

    def read_vectors_2d(self,plane):
        for wavenumberii in range(self.NumberWaveNumbers):
            buff=".%03d" % (wavenumberii+2)
            filename=PackagePath+"/qvectors/qvectors3d/qvector"+str(buff)
            qvector_file = Path(filename)
            if qvector_file.is_file():
                vector=pd.read_csv(filename,header=None,sep="\s+")
                if plane == "xy":
                    x_vector=vector.loc[:,0]
                    y_vector=vector.loc[:,1]
                    z_vector=pd.Series(np.zeros(len(x_vector)))
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
                elif plane == "yz":
                    y_vector=vector.loc[:,0]
                    x_vector=pd.Series(np.zeros(len(x_vector)))
                    z_vector=vector.loc[:,1]
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
                elif plane == "xz":
                    z_vector=vector.loc[:,1]
                    x_vector=vector.loc[:,0]
                    y_vector=pd.Series(np.zeros(len(x_vector)))
                    vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                    vector*=self.wavegrid_spacing
            else:
                print("ERROR:wave_vector::read_vectors_2d, "+filename+" can not open")
                self.wavevector.append(vector)    
    
    def read_vectors_3d(self):
        for wavenumberii in range(self.NumberWaveNumbers):
            buff=".%03d" % (wavenumberii+2)
            filename=PackagePath+"/qvectors/qvectors3d/qvector"+str(buff)
            qvector_file = Path(filename)
            if qvector_file.is_file():
                vector=pd.read_csv(filename,header=None,sep="\s+")
                x_vector=vector.loc[:,0]
                y_vector=vector.loc[:,1]
                z_vector=vector.loc[:,2]
                vector=pd.concat([x_vector,y_vector,z_vector],axis=1)
                vector*=self.wavegrid_spacing
                self.wavevector.append(vector)
            else:
                print("ERROR:wave_vector::read_vectors_3d, "+filename+" can not open")

    def calculate(self):
        pass

    def show_MeanWaveNumber(self,index):
        mean_wavenumber=0
        for vectorii in range(self.wavevector[index].shape[0]):
            mean_wavenumber+=np.linalg.norm(self.wavevector[index].loc[vectorii].values)
        mean_wavenumber/=self.wavevector[index].shape[0]
        return mean_wavenumber
    
    def show_StDevWaveNumber(self,index):
        mean_wavenumber=self.show_MeanWaveNumber(index)
        variance=0
        for vectorii in range(self.wavevector[index].shape[0]):
            variance+=np.power(np.linalg.norm(self.wavevector[index].loc[vectorii].values)-mean_wavenumber,2)
        variance/=self.wavevector[index].shape[0]
        stdev_=np.power(variance,0.5)
        return stdev_

    def show_MeanWaveVector(self,index):
        mean_wavevector=np.zeros(3)
        for vectorii in range(self.wavevector[index].shape[0]):
            mean_wavevector+=self.wavevector[index].loc[vectorii].values
        mean_wavevector/=self.wavevector[index].shape[0]
        return mean_wavevector

    def show_ApproxWaveNumber(self,index):
        return self.approx_wavenumber[index]

    def vectorlist(self,index):
        return self.wavevector[index]
    
    def vectorcount(self,index):
        return self.wavevector[index].shape[0]