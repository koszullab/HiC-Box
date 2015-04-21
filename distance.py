#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from numpy import diagonal, average
import numpy as np

def distance_diagonal_law(M):
    return np.array([average(diagonal(M,j)) for j in range(min(M.shape))])

def base_diagonal(M):
    L = distance_diagonal_law(M)
    return np.min(np.nonzero(L))

    
