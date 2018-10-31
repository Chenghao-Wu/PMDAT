#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
import pandas as pd

class MultiBody(object):
    def __init__(self,MultiBodyName,MultiBodyType,MultiBodyMethod,ListType,SpeciesName,*arg):
        self.MultiBodyName=MultiBodyName
        self.MultiBodyType=MultiBodyType
        self.MultiBodyMethod=MultiBodyMethod
        self.ListType=ListType
        self.SpeciesName=SpeciesName
        self.arg=arg
        self.c_data=0
        self.MultiBody_data=0
        self.c_vector=np.zeros(3)


class MultiBody_List(object):
    NumberMultiBody_count=0
    def __init__(self,System):
        self.System=System
        self.c_MultiBodies={}
        self.c_MultiBodyList={}
        self.c_data=0

    def create_MultiBodyList(self,MultiBodyName,MultiBodyType,MultiBodyMethod,ListType,SpeciesName,*arg):
        if len(arg[0]) % 2 == 0:
            self.c_MultiBodies[MultiBodyName]=MultiBody(MultiBodyName,MultiBodyType,MultiBodyMethod,ListType,SpeciesName,*arg)
            if MultiBodyMethod == "centroid":
                if ListType ==  "species_atomlist":
                    species_atomlist=arg[0][1::2]
                    species_typelist=arg[0][0::2]
                    if species_typelist==list(self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,SpeciesName,0,species_atomlist,:],:].loc[0,"type"]):
                        self.c_MultiBodies[MultiBodyName].c_data=self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,SpeciesName,:,species_atomlist,:],:].loc[:,["x","y","z"]]
                        MultiBodies_unwrap=self.System.show_unwrap(self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,SpeciesName,:,species_atomlist,:],:].copy())
                        self.c_MultiBodies[MultiBodyName].c_vector=MultiBodies_unwrap.loc[:,["x","y","z"]].groupby(level=[0,1,2]).nth(0)-MultiBodies_unwrap.loc[:,["x","y","z"]].groupby(level=[0,1,2]).nth(-1)
                        self.c_MultiBodies[MultiBodyName].c_vector.columns=["v_x","v_y","v_z"]
                    else:
                        print("ERROR:MultiBody::create_MultiBodyList, please set correct atom type")
                else:
                    print("ERROR:MultiBody::create_MultiBodyList, Please set suitable ListType")
                self.c_MultiBodies[MultiBodyName].MultiBody_data=self.calculate_CenterMassTrajectory(self.System.c_DataFrame.c_data.loc[pd.IndexSlice[:,SpeciesName,:,species_atomlist,:],:].copy())
            else:
                print("ERROR:MultiBody::create_MultiBodyList, Please set suitable MultiBodyMethod")
            MultiBody_List.NumberMultiBody_count+=1

        else:
            print("ERROR:MultiBody::create_MultiBodyList, Please set 'AtomType AtomIndex'")

    def combine_multibody_lists(self,*arg):
        self.c_MultiBodyList=arg
        if len(self.c_MultiBodyList)>1:
            MultiBodyList_CenterMassPos=[self.c_MultiBodies[MultiBodyii].MultiBody_data for MultiBodyii in self.c_MultiBodyList]
            MultiBodyList_Vector=[self.c_MultiBodies[MultiBodyii].c_vector for MultiBodyii in self.c_MultiBodyList]
        else:
            MultiBodyList_CenterMassPos=[self.c_MultiBodies[self.c_MultiBodyList[0]].MultiBody_data]
            MultiBodyList_Vector=[self.c_MultiBodies[self.c_MultiBodyList[0]].c_vector]
        MultiBodyFrame=pd.concat(MultiBodyList_CenterMassPos,axis=0)
        MultiBodyFrame.sort_index(level=0,axis=0,inplace=True)
        self.c_data=MultiBodyFrame.copy()
        MultiBodyFrame=pd.concat(MultiBodyList_Vector,axis=0)
        MultiBodyFrame.sort_index(level=0,axis=0,inplace=True)
        self.c_data=pd.concat([self.c_data,MultiBodyFrame],axis=1)

    def calculate_CenterMassTrajectory(self,dataframe):
        dataframe_unwrap=self.System.show_unwrap(dataframe)
        dataframe_unwrap.loc[:,"mass*x"]=dataframe_unwrap.loc[:,"mass"]*dataframe_unwrap.loc[:,"x"]
        dataframe_unwrap.loc[:,"mass*y"]=dataframe_unwrap.loc[:,"mass"]*dataframe_unwrap.loc[:,"y"]
        dataframe_unwrap.loc[:,"mass*z"]=dataframe_unwrap.loc[:,"mass"]*dataframe_unwrap.loc[:,"z"]
        dataframe_unwrap=dataframe_unwrap.loc[:,["x","y","z","mass","mass*x","mass*y","mass*z"]]
        Centroid_data=dataframe_unwrap.groupby(level=[0,1,2]).sum()
        
        Centroid_data.loc[:,"x"]=Centroid_data.loc[:,"mass*x"].values/Centroid_data.loc[:,"mass"].values
        Centroid_data.loc[:,"y"]=Centroid_data.loc[:,"mass*y"].values/Centroid_data.loc[:,"mass"].values
        Centroid_data.loc[:,"z"]=Centroid_data.loc[:,"mass*z"].values/Centroid_data.loc[:,"mass"].values
        Centroid_data.drop(["mass*x","mass*y","mass*z"],axis=1,inplace=True)
        return Centroid_data

    @property
    def unwrap_pos(self):
        self.check_if_combine_lists()
        c_data=self.c_data.loc[:,["x","y","z"]]
        c_data.astype("float64",copy=False)
        return c_data
    
    @property
    def vector(self):
        self.check_if_combine_lists()
        c_data=self.c_data.loc[:,["v_x","v_y","v_z"]]
        c_data.astype("float64",copy=False)
        return c_data

    @property
    def wrap_pos(self):
        self.check_if_combine_lists()
        c_data=self.c_data.loc[:,["x","y","z"]]
        c_data.astype("float64",copy=False)
        c_data.loc[:,"x"]=self.c_data.loc[:,"x"]-np.round(self.c_data.loc[:,"x"]/self.System.c_AbsBoxSize[0])*self.System.c_AbsBoxSize[0]
        c_data.loc[:,"y"]=self.c_data.loc[:,"y"]-np.round(self.c_data.loc[:,"y"]/self.System.c_AbsBoxSize[1])*self.System.c_AbsBoxSize[1]
        c_data.loc[:,"z"]=self.c_data.loc[:,"z"]-np.round(self.c_data.loc[:,"z"]/self.System.c_AbsBoxSize[2])*self.System.c_AbsBoxSize[2]
        return c_data

    def check_if_combine_lists(self):
        if isinstance(self.c_data,int):
            print("ERROR:MultiBody_List::unwrap_pos, please combine multibody list first!")
        else:
            pass