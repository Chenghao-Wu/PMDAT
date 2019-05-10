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

from SMMSAT.src.cython_func import create_list 
from SMMSAT.src.reader.Reader import *

class LAMMPSReader(Reader):

    def __init__(self,filename,start=0,mass={}):
        self.TrajectoryFormat="lammps"
        self.sys_NumberAtoms=0
        self.sys_NumberFrames=0
        self.start=start
        self.sys_data=pd.DataFrame([0])
        self.sys_BoxSize=[]
        self.sys_filename=filename
        self.set_BoxSize()
        self.read_NumberAtoms()
        self.read_NumberFrames()
        self.read_customformat()
        self.excute_Reader()
        self.sys_Mass=mass

    def set_BoxSize(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            with open(self.sys_filename) as f:
                reader=csv.reader(f)
                # index for statement
                index_start=-99
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    if line[0] == 'ITEM: BOX BOUNDS pp pp pp':
                        index_start=index
                    # read box size along x axis, parse the string to float
                    elif index == index_start+1:
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[1]))
                    # read box size along y axis, parse the string to float
                    elif index == index_start+2:
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[1]))
                    # read box size along z axis, parse the string to float
                    elif index == index_start+3:
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[0]))
                        self.sys_BoxSize.append(float(re.split('\s+',line[0])[1]))
                        return
                    else:
                        continue
        else:
            print("error,LAMMPSReader::set_BoxSize,file can not be found")


    def read_NumberAtoms(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            with open(self.sys_filename) as f:
                reader=csv.reader(f)
                index_start=-99
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    if line[0] == 'ITEM: NUMBER OF ATOMS':
                        index_start=index
                    # read box size along x axis, parse the string to float
                    elif index == index_start+1:
                        self.sys_NumberAtoms = int(line[0])
                        return None
        else:
            print("error,LAMMPSReader::read_NumberAtoms,file can not be found")

    def read_NumberFrames(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            f=open(self.sys_filename)
            n_lines = sum(1 for line in f)
            f.close()
            self.sys_NumberFrames=int(n_lines/(self.sys_NumberAtoms+9))
        else:
            print("error,LAMMPSReader::read_NumberFrames,file can not be found")

    def read_customformat(self):
        my_file = Path(self.sys_filename)
        if my_file.is_file():
            with open(self.sys_filename) as f:
                reader=csv.reader(f)
                # index for statement
                for index,line in enumerate(reader):
                    line[0]=line[0].strip()
                    if line[0][0:11] == 'ITEM: ATOMS':
                        self.sys_CustomFormat = line[0][11:].split()
                        return None

    def excute_Reader(self):
        colume_name=self.sys_CustomFormat
        skipframes=np.arange((self.sys_NumberAtoms+9)*self.start)
        skip_list=create_list.pandas_skiplist(np.array([0,1,2,3,4,5,6,7,8]),self.sys_NumberAtoms,self.sys_NumberFrames)
        skip_list=np.unique(list(skipframes)+list(skip_list))
        #print(skip_list)
        self.sys_NumberFrames=self.sys_NumberFrames-self.start
        data=pd.read_csv(self.sys_filename,header=None,sep="\s+",skiprows=list(skip_list))
        data=data.dropna(axis=1,how="any")
        data.columns=colume_name
        data[["x","y","z"]]=data[["x","y","z"]].astype("float32",copy=False)
        data[["type","ix","iy","iz"]]=data[["type","ix","iy","iz"]].astype("int32",copy=False)
        self.sys_data=data
        #default mass = 1
        #mass_list_frame=np.repeat(1.0,self.sys_NumberAtoms)
        #mass_list=np.tile(mass_list_frame,self.sys_NumberFrames)
        #self.sys_data.loc[:,"mass"]=mass_list
        #self.sys_data["mass"]=self.sys_data["mass"].astype("int32",copy=False)

    