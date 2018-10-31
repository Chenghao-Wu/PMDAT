#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import numpy as np
import pandas as pd
import csv
from pathlib import Path
import re

from cython_func import create_list 
from reader.Reader import *

class LAMMPSReader(Reader):

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
        my_file = Path(self.c_filename)
        if my_file.is_file():
            with open(self.c_filename) as f:
                reader=csv.reader(f)
                # index for statement
                index_start=-99
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    if line[0] == 'ITEM: BOX BOUNDS pp pp pp':
                        index_start=index
                    # read box size along x axis, parse the string to float
                    elif index == index_start+1:
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[1]))
                    # read box size along y axis, parse the string to float
                    elif index == index_start+2:
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[1]))
                    # read box size along z axis, parse the string to float
                    elif index == index_start+3:
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.c_BoxSize.append(float(re.split('\s+',line[0])[1]))
                        return
                    else:
                        continue
        else:
            print("error,LAMMPSReader::set_BoxSize,file can not be found")


    def read_NumberAtoms(self):
        my_file = Path(self.c_filename)
        if my_file.is_file():
            with open(self.c_filename) as f:
                reader=csv.reader(f)
                index_start=-99
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    if line[0] == 'ITEM: NUMBER OF ATOMS':
                        index_start=index
                    # read box size along x axis, parse the string to float
                    elif index == index_start+1:
                        self.c_NumberAtoms = int(line[0])
                        return None
        else:
            print("error,LAMMPSReader::read_NumberAtoms,file can not be found")

    def read_NumberFrames(self):
        my_file = Path(self.c_filename)
        if my_file.is_file():
            f=open(self.c_filename)
            n_lines = sum(1 for line in f)
            f.close()
            self.c_NumberFrames=int(n_lines/(self.c_NumberAtoms+9))
        else:
            print("error,LAMMPSReader::read_NumberFrames,file can not be found")

    def excute_Reader(self):
        colume_name=["type","x","y","z","ix","iy","iz"]
        skip_list=create_list.pandas_skiplist(np.array([0,1,2,3,4,5,6,7,8]),self.c_NumberAtoms,self.c_NumberFrames)
        data=pd.read_csv(self.c_filename,header=None,sep="\s+",skiprows=list(skip_list))
        data=data.dropna(axis=1,how="any")
        data.columns=colume_name
        data[["x","y","z"]].astype("float64",copy=False)
        self.c_data=data
        #default mass = 1
        mass_list_frame=np.repeat(1.0,self.c_NumberAtoms)
        mass_list=np.tile(mass_list_frame,self.c_NumberFrames)
        self.c_data.loc[:,"mass"]=mass_list

    