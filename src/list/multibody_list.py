#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
import pandas as pd
import os
from pathlib import Path
import gc
import sys

class MultiBody(object):
    def __init__(self,MultiBodyName,MultiBodyType,ListType,SpeciesName,*arg):
        self.MultiBodyName=MultiBodyName
        self.MultiBodyType=MultiBodyType
        self.ListType=ListType
        self.SpeciesName=SpeciesName
        self.atomlist=arg[0][1::2]
        self.typelist=arg[0][0::2]
        self.source_data=0
        self.sys_data=pd.DataFrame()
        self.sys_vector=np.zeros(3)
        
class MultiBodyList(object):
    NumberMultiBody_count=0
    def __init__(self,System):
        self.System=System
        self.ListType="MultibodyList"
        self.sys_MultiBodies={}
        self.sys_MultiBodyList={}
        self.sys_data=0

    def create_MultiBodyList(self,MultiBodyName,MultiBodyType,ListType,SpeciesName,*arg):
        self.check_MultibodyArg(arg)
        self.sys_MultiBodies[MultiBodyName]=MultiBody(MultiBodyName,MultiBodyType,ListType,SpeciesName,*arg)
        if ListType ==  "species_atomlist":
            species_atomlist=arg[0][1::2]
            species_typelist=arg[0][0::2]
            self.check_SpeciesType(SpeciesName,species_atomlist,species_typelist,arg)
            self.sys_MultiBodies[MultiBodyName].source_data=self.System.sys_data.loc[pd.IndexSlice[:,SpeciesName,:,species_atomlist,:],:].loc[:,["x","y","z","unwrap_x","unwrap_y","unwrap_z","mass"]]
            self.sys_MultiBodies[MultiBodyName].sys_data[["vector_x","vector_y","vector_z"]]=self.sys_MultiBodies[MultiBodyName].source_data.loc[:,["unwrap_x","unwrap_y","unwrap_z"]].groupby(level=[0,1,2]).nth(0)-self.sys_MultiBodies[MultiBodyName].source_data.loc[:,["unwrap_x","unwrap_y","unwrap_z"]].groupby(level=[0,1,2]).nth(-1)
            self.sys_MultiBodies[MultiBodyName].source_data[["mass*unwrap_x","mass*unwrap_y","mass*unwrap_z"]]=self.sys_MultiBodies[MultiBodyName].source_data[["unwrap_x","unwrap_y","unwrap_z"]].multiply(self.sys_MultiBodies[MultiBodyName].source_data["mass"],axis="index")
            self.sys_MultiBodies[MultiBodyName].sys_data[["unwrap_x","unwrap_y","unwrap_z","mass"]]=self.sys_MultiBodies[MultiBodyName].source_data[["mass*unwrap_x","mass*unwrap_y","mass*unwrap_z","mass"]].groupby(level=[0,1,2]).sum().divide(self.sys_MultiBodies[MultiBodyName].source_data["mass"].groupby(level=[0,1,2]).sum(),axis="index")
            self.sys_MultiBodies[MultiBodyName].sys_data[["ix","iy","iz"]]=np.round((self.sys_MultiBodies[MultiBodyName].sys_data[["unwrap_x","unwrap_y","unwrap_z"]]-(self.System.get_BoxSize[[1,3,5]]-self.System.get_AbsBoxSize[[0,1,2]]/2))/self.System.get_AbsBoxSize[[0,1,2]])
            ModificationSymmetricAxis=np.array([(self.System.get_BoxSize[1]+self.System.get_BoxSize[0])/2,(self.System.get_BoxSize[3]+self.System.get_BoxSize[2])/2,(self.System.get_BoxSize[5]+self.System.get_BoxSize[4])/2])
            self.sys_MultiBodies[MultiBodyName].sys_data[["x","y","z"]]=pd.DataFrame(self.sys_MultiBodies[MultiBodyName].sys_data[["unwrap_x","unwrap_y","unwrap_z"]].values-self.sys_MultiBodies[MultiBodyName].sys_data[["ix","iy","iz"]].values*self.System.get_AbsBoxSize[[0,1,2]]+ModificationSymmetricAxis,columns=["x","y","z"],index=self.sys_MultiBodies[MultiBodyName].sys_data.index)
        else:
            print("ERROR:MultiBody::create_MultiBodyList, Please set suitable ListType")
        MultiBodyList.NumberMultiBody_count+=1

    #index: frame, speciesname, speciesindex,(slice any multibody:.loc[frameindex].loc[speciesname].loc[speciesindex].iloc[multibody index])
    def combine_multibody_lists(self,*arg):
        self.sys_MultiBodyList=arg
        if isinstance(arg,tuple):
            MultiBodyList_sys_data=[self.sys_MultiBodies[multibodyii].sys_data for multibodyii in arg]
            self.sys_data=pd.concat(MultiBodyList_sys_data,axis=0)
            self.sys_data.sort_index(level=0,axis=0,inplace=True)
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]
        else:
            multibody=self.sys_MultiBodyList[0]
            self.sys_data=self.sys_MultiBodies[multibody].sys_data
            self.wrap_pos=self.sys_data[["x","y","z"]]
            self.unwrap_pos=self.sys_data[["unwrap_x","unwrap_y","unwrap_z"]]

    #check setting
    def check_MultibodyArg(self,arg):
        if len(arg[0]) % 2 == 0:
            pass
        else:
            sys.exit(2)

    def check_SpeciesType(self,SpeciesName,species_atomlist,species_typelist,arg):
        
        if species_typelist==list(self.System.sys_data.loc[pd.IndexSlice[:,SpeciesName,0,species_atomlist,:],:].loc[0,"type"]):
            pass
        else:
            sys.exit(2)

    def check_MultibodySetting(self,Multibody):
        if Multibody.typelist==list(self.System.sys_data.loc[pd.IndexSlice[:,Multibody.SpeciesName,0,Multibody.atomlist,:],:].loc[0,"type"]):
            return True
        else:
            return False

    def check_if_combine_lists(self):
        if isinstance(self.sys_data,int):
            print("ERROR:MultiBody_List::unwrap_pos, please combine multibody list first!")
        else:
            pass
