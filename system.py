#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""

import pandas as pd
import numpy as np
from analysis import *

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
    def __init__(self,DataFrame):
        self.c_NumberBlocks     =   int(0)
        self.c_BlockSize        =   int(0)
        self.SpeciesList        =   {}
        self.c_DataFrame        =   DataFrame
        self.c_NumberSteps      =   int(0)
        self.c_ExponentialBase  =   int(0)
        self.c_TimeUnit         =   float(0.0)
        self.c_TimeList         =   np.array([0])
        self.c_NumberTimeGaps   =   0
        self.c_TimeGap          =   np.zeros(1)
        self.c_AbsBoxSize       =   [self.c_DataFrame.c_BoxSize[1]-self.c_DataFrame.c_BoxSize[0],self.c_DataFrame.c_BoxSize[3]-self.c_DataFrame.c_BoxSize[2],self.c_DataFrame.c_BoxSize[5]-self.c_DataFrame.c_BoxSize[4]]

    def get_NumberBlocks(self):
        return self.c_NumberBlocks
    
    def get_Data(self):
        return self.c_DataFrame.c_data

    def get_NumberAtoms(self):
        return self.c_DataFrame.c_NumberAtoms

    def get_NumberFrames(self):
        return self.c_DataFrame.c_NumberFrames
    
    def get_BoxSize(self):
        return self.c_DataFrame.c_BoxSize

    def show_n_timesteps(self):
        return self.c_NumberSteps

    def set_LinearTimeScheme(self,n_blocks,TimeUnit):
        self.c_NumberSteps=n_blocks
        self.c_NumberBlocks=n_blocks-1
        self.c_TimeUnit=TimeUnit
        self.c_BlockSize=1
        self.c_NumberTimeGaps   =   self.c_BlockSize + self.c_NumberBlocks

    def set_ExponentialTimeScheme(self,n_blocks,blocksize,ExponentialBase,TimeUnit):
        self.c_NumberBlocks=n_blocks
        self.c_ExponentialBase=ExponentialBase
        self.c_TimeUnit=TimeUnit
        self.c_BlockSize=blocksize
        self.c_NumberSteps=self.c_NumberBlocks*self.c_BlockSize+1
        self.c_NumberTimeGaps   =   self.c_BlockSize + self.c_NumberBlocks

    
    def set_Species(self,SpeciesName,NumberSpecies,Atomtype,Atoms):# molecule setting
        self.SpeciesList[SpeciesName]=Species(SpeciesName,NumberSpecies,Atomtype,Atoms)
        System.Species_count+=1

    def check_SpeciesSetting(self):
        total_atoms=0
        atom_type=[]
        for speciessii in self.SpeciesList:
            total_atoms+=sum(self.SpeciesList[speciessii].AtomsList)*self.SpeciesList[speciessii].NumberSpecies
            atom_type=self.SpeciesList[speciessii].Atomtype
        if total_atoms == self.c_DataFrame.c_NumberAtoms:
            print("\nINIT: SMMAT Read "+ str(total_atoms)+ " Atoms")
        else:
            print("\nERROR:System::check_SpeciesSetting, Numer of Atoms in Trj file is not the same with SpeciesSetting")
        for atomtypeii in atom_type:
            if atomtypeii in self.c_DataFrame.c_data.loc[:,"type"].values:
                continue
            else:
                print("\nERROR:System::check_SpeciesSetting, "+atomtypeii+" is not in Trj file")
                break
    
    def create_TimeList(self):
        timeii=0
        block_starttime=0
        self.c_TimeList=np.zeros(self.c_NumberSteps,dtype=np.float64)
        for blockii in range(0,self.c_NumberBlocks):
            for expstepii in range(1,self.c_BlockSize+1):
                timeii+=1
                if np.power(self.c_ExponentialBase,expstepii-1) <= expstepii:
                    self.c_TimeList[timeii] = block_starttime+float(expstepii)*self.c_TimeUnit
                else:
                    self.c_TimeList[timeii] = block_starttime+float(np.floor(np.power(self.c_ExponentialBase,expstepii-1)))*self.c_TimeUnit
            block_starttime = self.c_TimeList[timeii]
    
    def check_TimeScheme(self):
        if self.c_NumberSteps == self.c_DataFrame.c_NumberFrames:
            print("\nINIT: SMMAT Read "+str(self.c_NumberSteps)+" Frames")
        else:
            print("\nERROR:System::check_TimeScheme, timescheme setting is not correct")

    def set_DataFrameIndex(self):

        AtomIndex_frame=np.arange(self.c_DataFrame.c_NumberAtoms)
        AtomIndex=np.tile(AtomIndex_frame,self.c_DataFrame.c_NumberFrames)

        SpeciesTypeIndex_frame=[]
        for SpeciesTypeii in self.SpeciesList:
            SpeciesTypeLength=np.sum(self.SpeciesList[SpeciesTypeii].AtomsList)*self.SpeciesList[SpeciesTypeii].NumberSpecies
            SpeciesTypeIndex_frame.append(np.repeat([self.SpeciesList[SpeciesTypeii].SpeciesName],SpeciesTypeLength))

        SpeciesTypeIndex_frame=np.concatenate(SpeciesTypeIndex_frame,axis=0)
        
        SpeciesTypeIndex=np.tile(SpeciesTypeIndex_frame,self.c_DataFrame.c_NumberFrames)

        SpeciesIndex_frame=[]
        for SpeciesTypeii in self.SpeciesList:
            for Speciesii in range(self.SpeciesList[SpeciesTypeii].NumberSpecies):
                SpeciesLength=np.sum(self.SpeciesList[SpeciesTypeii].AtomsList)
                SpeciesIndex_frame.append(np.repeat([Speciesii],SpeciesLength))
            #print(SpeciesIndex_frame,len(SpeciesIndex_frame))
        
        SpeciesIndex_frame=np.concatenate(SpeciesIndex_frame,axis=0)
        
        SpeciesIndex=np.tile(SpeciesIndex_frame,self.c_DataFrame.c_NumberFrames)
        
        SpeciesAtomIndex_frame=[]
        for SpeciesTypeii in self.SpeciesList:
            for Speciesii in range(self.SpeciesList[SpeciesTypeii].NumberSpecies):
                SpeciesAtomIndex_frame.append(np.arange(np.sum(self.SpeciesList[SpeciesTypeii].AtomsList)))
        SpeciesAtomIndex_frame=np.concatenate(SpeciesAtomIndex_frame,axis=0)
        SpeciesAtomIndex=np.tile(SpeciesAtomIndex_frame,self.c_DataFrame.c_NumberFrames)


        FrameIndex=np.repeat(np.arange(0,self.c_DataFrame.c_NumberFrames),self.c_DataFrame.c_NumberAtoms)
        Index_5d=pd.MultiIndex.from_arrays([FrameIndex,SpeciesTypeIndex,SpeciesIndex,SpeciesAtomIndex,AtomIndex])
        self.c_DataFrame.c_data.index=Index_5d

        
    def calculate_TimeGap(self):
        self.c_TimeGap=np.zeros( self.c_NumberTimeGaps, dtype=np.float64)
        for timeii in range(self.c_BlockSize):
            self.c_TimeGap[timeii]=self.c_TimeList[timeii]-self.c_TimeList[0]
        for timeii in range(1,self.c_NumberBlocks):
            self.c_TimeGap[timeii+self.c_BlockSize-1]=self.c_TimeList[self.c_BlockSize*timeii]-self.c_TimeList[0]
        self.c_TimeGap[self.c_NumberTimeGaps-1]=self.c_TimeList[self.c_NumberSteps-1]-self.c_TimeList[0]

    def TimeGapWeighting(self,fullblock):
        weighting=np.zeros(self.c_NumberTimeGaps)
        for timegapii in range(self.c_BlockSize):
            weighting[timegapii]=self.c_NumberBlocks
        for timegapii in range(self.c_BlockSize,self.c_NumberTimeGaps-1):
            block_timegapii = timegapii - self.c_BlockSize+1
            weighting[timegapii]=(self.c_NumberBlocks-block_timegapii)+int(fullblock)*(self.c_NumberBlocks-block_timegapii)*(self.c_BlockSize-1)+1
        weighting[self.c_NumberTimeGaps-1]=1
        return weighting

    def convert_TimeToFrame(self,Time):
        FrameIndex=np.int(Time/self.c_TimeUnit)
        return FrameIndex

    def show_NumberFrame(self,TimeGapii):
        FrameIndex=self.convert_TimeToFrame(TimeGapii)

    def show_unwrap(self,dataframe):
        dataframe.loc[:,"x"]=dataframe.loc[:,"x"]+np.sum(np.absolute(self.c_DataFrame.c_BoxSize[1]-self.c_DataFrame.c_BoxSize[0]))*dataframe.loc[:,"ix"]
        dataframe.loc[:,"y"]=dataframe.loc[:,"y"]+np.sum(np.absolute(self.c_DataFrame.c_BoxSize[3]-self.c_DataFrame.c_BoxSize[2]))*dataframe.loc[:,"iy"]
        dataframe.loc[:,"z"]=dataframe.loc[:,"z"]+np.sum(np.absolute(self.c_DataFrame.c_BoxSize[5]-self.c_DataFrame.c_BoxSize[4]))*dataframe.loc[:,"iz"]
        return dataframe

    def show_wrap(self,dataframe):
        return dataframe

    def convert_WrapToUnwrap(self,dataframe):
        pass
