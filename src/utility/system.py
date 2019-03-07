#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import pandas as pd
import numpy as np


class Species(object):
    def __init__(self,SpeciesName,NumberSpecies,Atomtype,Atoms):
        self.SpeciesName=SpeciesName
        self.NumberSpecies=NumberSpecies
        self.Atomtype=Atomtype
        self.AtomsList=Atoms
        

class System(object):
    MAX_NumberSpecies=10
    displacement_limit=0
    Species_count=0
    Version="v0.1"


#   system parameters
    def __init__(self,Reader):
        self.sys_NumberBlocks     =   int(0)
        self.sys_BlockSize        =   int(0)
        self.sys_SpeciesDict      =   {}
        self.sys_Reader           =   Reader
        self.sys_NumberSteps      =   int(0)
        self.sys_ExponentialBase  =   int(0)
        self.sys_TimeUnit         =   float(0.0)
        self.sys_TimeList         =   np.array([0])
        self.sys_NumberTimeGaps   =   0
        self.sys_TimeSchemeType   =   0
        self.sys_TimeGapTable          =   np.zeros(1)
        self.sys_AbsBoxSize       =   np.array([self.sys_Reader.sys_BoxSize[1]-self.sys_Reader.sys_BoxSize[0],self.sys_Reader.sys_BoxSize[3]-self.sys_Reader.sys_BoxSize[2],self.sys_Reader.sys_BoxSize[5]-self.sys_Reader.sys_BoxSize[4]])
        self.sys_data=Reader.sys_data
        self.build_UnwrapPos()
        

#   get function for properties in analysis systems
    @property
    def get_NumberBlocks(self):
        return self.sys_NumberBlocks

    @property
    def get_BlockSize(self):
        return self.sys_BlockSize

    @property
    def get_Data(self):
        return self.sys_Reader.sys_data

    @property
    def get_NumberAtoms(self):
        return self.sys_Reader.sys_NumberAtoms

    @property
    def get_NumberFrames(self):
        return self.sys_Reader.sys_NumberFrames
    
    @property
    def get_BoxSize(self):
        return np.array(self.sys_Reader.sys_BoxSize)

    @property
    def get_AbsBoxSize(self):
        return self.sys_AbsBoxSize

    @property
    def get_NumberTimesteps(self):
        return self.sys_NumberSteps
    
    @property
    def get_NumberTimeGaps(self):
        return self.sys_NumberTimeGaps
    
    @property
    def get_TimeGapTable(self):
        return self.sys_TimeGapTable

    @property
    def get_NumberSteps(self):
        return self.sys_NumberSteps

    @property
    def get_TimeList(self):
        return self.sys_TimeList

#   set function to obtain parameters from user
    def set_LinearTimeScheme(self,n_blocks,TimeUnit):
        self.sys_NumberSteps=n_blocks
        self.sys_NumberBlocks=n_blocks-1
        self.sys_TimeUnit=TimeUnit
        self.sys_BlockSize=1
        self.sys_NumberTimeGaps   =    self.sys_NumberBlocks
        self.sys_TimeSchemeType="linear"

    def set_ExponentialTimeScheme(self,n_blocks,blocksize,ExponentialBase,TimeUnit):
        self.sys_NumberBlocks=n_blocks
        self.sys_ExponentialBase=ExponentialBase
        self.sys_TimeUnit=TimeUnit
        self.sys_BlockSize=blocksize
        self.sys_NumberSteps=self.sys_NumberBlocks*self.sys_BlockSize+1
        self.sys_NumberTimeGaps   =   self.sys_BlockSize + self.sys_NumberBlocks
        self.sys_TimeSchemeType="exponential"
    
    def set_Species(self,SpeciesName,NumberSpecies,Atomtype,Atoms):# molecule setting
        self.sys_SpeciesDict[SpeciesName]=Species(SpeciesName,NumberSpecies,Atomtype,Atoms)
        System.Species_count+=1

#   check function check input parameters from users
    def check_TimeScheme(self):
        if self.sys_NumberSteps == self.sys_Reader.sys_NumberFrames:
            print("\nINIT: SMMAT Read "+str(self.sys_NumberSteps)+" Frames")
        else:
            print("\nERROR:System::check_TimeScheme, timescheme setting is not correct")

    def check_SpeciesSetting(self):
        total_atoms=0
        atom_type=[]
        for speciessii in self.sys_SpeciesDict:
            total_atoms+=sum(self.sys_SpeciesDict[speciessii].AtomsList)*self.sys_SpeciesDict[speciessii].NumberSpecies
            atom_type=self.sys_SpeciesDict[speciessii].Atomtype
        if total_atoms == self.sys_Reader.sys_NumberAtoms:
            print("\nINIT: SMMAT Read "+ str(total_atoms)+ " Atoms")
        else:
            print("\nERROR:System::check_SpeciesSetting, Numer of Atoms in Trj file is not the same with SpeciesSetting")
        for atomtypeii in atom_type:
            if atomtypeii in self.sys_Reader.sys_data.loc[:,"type"].values:
                continue
            else:
                print("\nERROR:System::check_SpeciesSetting, "+atomtypeii+" is not in Trj file")
                break

    def check_TrjFormat(self):
        columns=self.sys_Reader.sys_data.columns
        if "ix" in columns and "iy" in columns and "iz" in columns:
            pass
        else:
            image=np.repeat(np.repeat(0,self.sys_Reader.sys_NumberAtoms),self.sys_Reader.sys_NumberFrames)
            self.sys_Reader.sys_data.loc[:,["ix","iy","iz"]]=image

    def build_UnwrapPos(self):
        self.sys_data["unwrap_x"]=self.sys_data["ix"]*self.sys_AbsBoxSize[0]+self.sys_data["x"]
        self.sys_data["unwrap_y"]=self.sys_data["iy"]*self.sys_AbsBoxSize[1]+self.sys_data["y"]
        self.sys_data["unwrap_z"]=self.sys_data["iz"]*self.sys_AbsBoxSize[2]+self.sys_data["z"]
        
#   build function to build lists
    def build_TimeList(self):
        timeii=0
        block_starttime=0
        self.sys_TimeList=np.zeros(self.sys_NumberSteps,dtype=np.float32)
        for blockii in range(0,self.sys_NumberBlocks):
            for expstepii in range(1,self.sys_BlockSize+1):
                timeii+=1
                if np.power(self.sys_ExponentialBase,expstepii-1) <= expstepii:
                    self.sys_TimeList[timeii] = block_starttime+float(expstepii)*self.sys_TimeUnit
                else:
                    self.sys_TimeList[timeii] = block_starttime+float(np.floor(np.power(self.sys_ExponentialBase,expstepii-1)))*self.sys_TimeUnit
            block_starttime = self.sys_TimeList[timeii]

    def build_DataFrameIndex(self):
        AtomIndex_frame=np.arange(self.sys_Reader.sys_NumberAtoms)
        AtomIndex=np.tile(AtomIndex_frame,self.sys_Reader.sys_NumberFrames)
        AtomIndex=AtomIndex.astype("int32")
        SpeciesTypeIndex_frame=[]
        for SpeciesTypeii in self.sys_SpeciesDict:
            SpeciesTypeLength=np.sum(self.sys_SpeciesDict[SpeciesTypeii].AtomsList)*self.sys_SpeciesDict[SpeciesTypeii].NumberSpecies
            SpeciesTypeIndex_frame.append(np.repeat([self.sys_SpeciesDict[SpeciesTypeii].SpeciesName],SpeciesTypeLength))

        SpeciesTypeIndex_frame=np.concatenate(SpeciesTypeIndex_frame,axis=0)
        
        SpeciesTypeIndex=np.tile(SpeciesTypeIndex_frame,self.sys_Reader.sys_NumberFrames)

        SpeciesIndex_frame=[]
        for SpeciesTypeii in self.sys_SpeciesDict:
            for Speciesii in range(self.sys_SpeciesDict[SpeciesTypeii].NumberSpecies):
                SpeciesLength=np.sum(self.sys_SpeciesDict[SpeciesTypeii].AtomsList)
                SpeciesIndex_frame.append(np.repeat([Speciesii],SpeciesLength))
        
        SpeciesIndex_frame=np.concatenate(SpeciesIndex_frame,axis=0)
        SpeciesIndex_frame=SpeciesIndex_frame.astype("int32")
        
        SpeciesIndex=np.tile(SpeciesIndex_frame,self.sys_Reader.sys_NumberFrames)
        SpeciesAtomIndex_frame=[]
        for SpeciesTypeii in self.sys_SpeciesDict:
            for Speciesii in range(self.sys_SpeciesDict[SpeciesTypeii].NumberSpecies):
                SpeciesAtomIndex_frame.append(np.arange(np.sum(self.sys_SpeciesDict[SpeciesTypeii].AtomsList)))
        SpeciesAtomIndex_frame=np.concatenate(SpeciesAtomIndex_frame,axis=0)
        SpeciesAtomIndex=np.tile(SpeciesAtomIndex_frame,self.sys_Reader.sys_NumberFrames)
        SpeciesAtomIndex=SpeciesAtomIndex.astype("int32")

        FrameIndex=np.repeat(np.arange(0,self.sys_Reader.sys_NumberFrames),self.sys_Reader.sys_NumberAtoms)
        FrameIndex=FrameIndex.astype("int32")
        Index_5d=pd.MultiIndex.from_arrays([FrameIndex,SpeciesTypeIndex,SpeciesIndex,SpeciesAtomIndex,AtomIndex])
        self.sys_data.index=Index_5d
        
    def build_TimeGapList(self):
        if self.sys_NumberTimeGaps>0:
            self.sys_TimeGapTable=np.zeros( self.sys_NumberTimeGaps, dtype=np.float32)
            for timeii in range(self.sys_BlockSize):
                self.sys_TimeGapTable[timeii]=self.sys_TimeList[timeii]-self.sys_TimeList[0]
            for timeii in range(1,self.sys_NumberBlocks):
                self.sys_TimeGapTable[timeii+self.sys_BlockSize-1]=self.sys_TimeList[self.sys_BlockSize*timeii]-self.sys_TimeList[0]
                self.sys_TimeGapTable[self.sys_NumberTimeGaps-1]=self.sys_TimeList[self.sys_NumberSteps-1]-self.sys_TimeList[0]

    def convert_TimeToFrame(self,Time):
        FrameIndex=np.int(Time/self.sys_TimeUnit)
        return FrameIndex

    def show_NumberFrame(self,TimeGapii):
        FrameIndex=self.convert_TimeToFrame(TimeGapii)

    def show_unwrap(self,dataframe):
        dataframe.loc[:,"x"]=dataframe.loc[:,"x"]+np.sum(np.absolute(self.sys_Reader.sys_BoxSize[1]-self.sys_Reader.sys_BoxSize[0]))*dataframe.loc[:,"ix"]
        dataframe.loc[:,"y"]=dataframe.loc[:,"y"]+np.sum(np.absolute(self.sys_Reader.sys_BoxSize[3]-self.sys_Reader.sys_BoxSize[2]))*dataframe.loc[:,"iy"]
        dataframe.loc[:,"z"]=dataframe.loc[:,"z"]+np.sum(np.absolute(self.sys_Reader.sys_BoxSize[5]-self.sys_Reader.sys_BoxSize[4]))*dataframe.loc[:,"iz"]
        return dataframe

    def show_wrap(self,dataframe):
        return dataframe

    def convert_WrapToUnwrap(self,dataframe):
        pass
    
    



