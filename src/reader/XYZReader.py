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

from SMMSAT.src.reader import create_list 
from SMMSAT.src.reader.Reader import *

class XYZReader(Reader):

    def __init__(self,filename,LogFilename=None):
        self.c_NumberAtoms=0
        self.c_NumberFrames=0
        self.c_data=pd.DataFrame([0])
        self.c_BoxSize=[]
        self.c_filename=filename
        self.set_BoxSize()
        self.read_NumberAtoms()
        self.read_NumberFrames()
        self.excute_Reader()
        self.LogFilename=LogFilename
        self.ReadLog=False
        
    def set_BoxSize(self):
        pass

    def read_NumberAtoms(self):
        my_file = Path(self.c_filename)
        if my_file.is_file():
            with open(self.c_filename) as f:
                reader=csv.reader(f)
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    self.c_NumberAtoms=int(line[0])
                    return None
        else:
            print("error,XYZReader::read_NumberAtoms,file can not be found")

    def read_NumberFrames(self):
        my_file = Path(self.c_filename)
        if my_file.is_file():
            f=open(self.c_filename)
            n_lines = sum(1 for line in f)
            f.close()
            self.c_NumberFrames=int(n_lines/(self.c_NumberAtoms+2))
        else:
            print("error,XYZReader::read_NumberFrames,file can not be found")

    def excute_Reader(self):
        colume_name=["type","x","y","z"]
        skip_list=create_list.pandas_skiplist(np.array([0,1]),self.c_NumberAtoms,self.c_NumberFrames)
        data=pd.read_csv(self.c_filename,header=None,sep="\s+",skiprows=skip_list)
        data=data.dropna(axis=1,how="any")
        data.columns=colume_name
        data[["x","y","z"]].astype("float64",copy=False)
        self.c_data=data
        mass_list_frame=np.repeat(1.0,self.c_NumberAtoms)
        mass_list=np.tile(mass_list_frame,self.c_NumberFrames)
        self.c_data.loc[:,"mass"]=mass_list