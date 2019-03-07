#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:54:59 2018

@author: bruce
"""
import numpy as np
from scipy.special import legendre

data_precision=6

def unit_vector(vector):
    if vector.shape == [3,]:
        return vector / np.linalg.norm(vector)
    else:
        return vector / np.linalg.norm(vector,axis=1).reshape(len(vector),1)

def legendre_polynomial(array,level):
    return np.polyval(legendre(level),array)