#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import numpy as np
import pandas as pd

class CreateList(object):
    def __init__(self,System):
        self.System=System
        
    def create_list(self,*keyword):
        if keyword[0]=="all":
            self.c_data=self.System.c_DataFrame.c_data
        elif keyword[0]=="type_system":
            #Selects atoms of type atom_type
            self.AtomType=keyword[1]
            self.c_data=self.System.c_DataFrame.c_data.loc[self.System.c_DataFrame.c_data["type"]==self.AtomType]

        elif keyword[0]=="type_species":
            #Selects atoms of type atom_type in species species_name
            self.SpeciesName=keyword[1]
            self.AtomType=keyword[2]
            self.c_data=self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,self.SpeciesName,:,:,:],:].loc[self.System.c_DataFrame.c_data["type"]==self.AtomType]
        elif keyword[0]=="species":
            #Selects all atoms within species species_name
            self.SpeciesName=keyword[1]
            self.c_data=self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,self.SpeciesName,:,:,:],:]
        else:
            print("\nERROR:List::create_list, Please set correct keyword\n")
    
    @property
    def unwrap_pos(self):
        c_data=self.c_data.loc[:,["x","y","z"]]
        c_data.astype("float64",copy=False)
        c_data.loc[:,"x"]=self.c_data.loc[:,"x"]+self.System.c_AbsBoxSize[0]*self.System.c_DataFrame.c_data.loc[:,"ix"]
        c_data.loc[:,"y"]=self.c_data.loc[:,"y"]+self.System.c_AbsBoxSize[1]*self.System.c_DataFrame.c_data.loc[:,"iy"]
        c_data.loc[:,"z"]=self.c_data.loc[:,"z"]+self.System.c_AbsBoxSize[2]*self.System.c_DataFrame.c_data.loc[:,"iz"]
        return c_data
    
    @property
    def wrap_pos(self):
        c_data=self.c_data.loc[:,["x","y","z"]]
        c_data.astype("float64",copy=False)
        return c_data
