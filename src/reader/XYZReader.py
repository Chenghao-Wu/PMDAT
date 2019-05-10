#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
import pandas as pd
import csv
import pathlib

from SMMSAT.src.cython_func import create_list 
from SMMSAT.src.reader.Reader import *

class XYZReader(Reader):

    def __init__(self,filename,box=[0,0,0,0,0,0],mass={},Type=None): #box[xlo,xhi,ylo,yhi,zlo,zhi]
        self.TrajectoryFormat="xyz"
        self.sys_NumberAtoms=0
        self.sys_NumberFrames=0
        self.sys_data=pd.DataFrame([0])
        self.sys_filename=filename
        self.set_BoxSize()
        self.read_NumberAtoms()
        self.read_NumberFrames()
        self.excute_Reader()
        self.sys_BoxSize=np.array(box)
        self.mass=mass
        self.ReadLog=False
        self.Type=Type

    def read_NumberAtoms(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            with open(self.sys_filename) as f:
                reader=csv.reader(f)
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    self.sys_NumberAtoms=int(line[0])
                    return None
        else:
            print("error,XYZReader::read_NumberAtoms,file can not be found")

    def read_NumberFrames(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            f=open(self.sys_filename)
            n_lines = sum(1 for line in f)
            f.close()
            self.sys_NumberFrames=int(n_lines/(self.sys_NumberAtoms+2))
        else:
            print("error,XYZReader::read_NumberFrames,file can not be found")

    def excute_Reader(self):
        colume_name=["type","x","y","z"]
        skip_list=create_list.pandas_skiplist(np.array([0,1]),self.sys_NumberAtoms,self.sys_NumberFrames)
        data=pd.read_csv(self.sys_filename,header=None,sep="\s+",skiprows=skip_list)
        data=data.dropna(axis=1,how="any")
        data.columns=colume_name
        data[["x","y","z"]].astype("float64",copy=False)
        image=np.repeat(0,self.sys_NumberAtoms)
        image_list=np.tile(image,self.sys_NumberFrames)
        data["ix"]=image_list
        data["iy"]=image_list
        data["iz"]=image_list
        #data["type"]=1
        self.sys_data=data
        mass_list_frame=np.repeat(1.0,self.sys_NumberAtoms)
        mass_list=np.tile(mass_list_frame,self.sys_NumberFrames)
        self.sys_data.loc[:,"mass"]=mass_list
        
        