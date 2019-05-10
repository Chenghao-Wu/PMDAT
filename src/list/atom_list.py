#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import numpy as np
import pandas as pd

class AtomList(object):
    def __init__(self,System,*keyword):
        self.System=System
        self.ListType="AtomList"
        if keyword[0]=="all":
            self.AtomListType="all"
            self.sys_data=self.System.sys_data
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
        elif keyword[0]=="type_system":
#           Selects atoms of type atom_type
            self.AtomListType="type_system"
            self.AtomType=keyword[1]
            self.sys_data=self.System.sys_data.loc[self.System.sys_data["type"]==self.AtomType]
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
        elif keyword[0]=="type_species":
#           Selects atoms of type atom_type in species species_name
            self.AtomListType="type_species"
            self.SpeciesName=keyword[1]
            if isinstance(keyword[2],list):
                self.AtomType=keyword[2]
            else:
                self.AtomType=[keyword[2]]
            atomselection=0
            for type_speciessii in self.AtomType:
                atomselection=np.logical_or(self.System.sys_data["type"]==type_speciessii,atomselection)
            self.sys_data=self.System.sys_data.loc[pd.IndexSlice[:,self.SpeciesName,:,:,:],:].loc[atomselection]
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
            self.selectedspecieslength=self.sys_data.loc[pd.IndexSlice[0,self.SpeciesName,0,:,:],:].shape[0]
        elif keyword[0]=="species":
#           Selects all atoms within species species_name
            self.AtomListType="species"
            self.SpeciesName=keyword[1]
            self.sys_data=self.System.sys_data.loc[pd.IndexSlice[:,self.SpeciesName,:,:,:],:]
            self.sys_NumberAtomPerChain=sum(self.System.sys_SpeciesDict[self.SpeciesName].AtomsList)
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
            self.selectedspecieslength=self.sys_data.loc[pd.IndexSlice[0,self.SpeciesName,0,:,:],:].shape[0]
        elif keyword[0] == "index_atom":
            self.AtomListType="index_atom"
            self.SpeciesName=keyword[1]
            self.index_atom=keyword[2]
            self.sys_data=self.System.sys_data.loc[pd.IndexSlice[:,self.SpeciesName,:,self.index_atom,:],:]
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
        else:
            print("\nERROR:List::create_list, Please set correct keyword\n")
        