#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import xml.etree.ElementTree as ET 
from pathlib import Path
import numpy as np

class Reader(object):
    def __init__(self,filename):
        pass

    def excute(self):
        pass

    def excute_Reader(self):
        pass

    def set_BoxSize(self):
        pass
    def read_NumberAtoms(self):
        pass
    def read_NumberFrames(self):
        pass
    
    def read_xml(self,filename,arg):
        tree = ET.parse(filename)
        root = tree.getroot()
        arg_list=[argii.text for argii in root.iter(arg)][0].split("\n")[1:-1]
        return arg_list
    """
    def set_ReadLog(self,*arg):
        
        Log file is xml format
        
        if self.LogFilename==None:
            print("ERROR:Reader::set_ReadLog, Please set the LogFile Name")
            self.ReadLog=True
        for argii in arg:
            if argii=="all":
                if len(arg)==1:
                    pass
                else:
                    print("Error:Reader::set_ReadLog, other parameters are not allowed to be used together with 'all'")
            if argii=="mass":
                mass_list_frame=np.array(self.read_xml(self.LogFilename,"mass"),dtype=np.float64)
                mass_list=np.tile(mass_list_frame,self.c_NumberFrames)
                self.c_data.loc[:,"mass"]=mass_list
            if argii == "box":
                self.c_BoxSize=np.array(self.read_xml(self.LogFilename,"box")[0].split(),dtype=np.float64)
    """