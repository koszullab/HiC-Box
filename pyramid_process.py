#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from mirnylib import numutils as ntls
from numpy import log,log10,log2
from numpy.ma import corrcoef
from scipy import ndimage
import distance
import numpy as np
import matplotlib.pyplot as plt

cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))

def process(t,matrix_name,normalization,order,iterations,exposant,gaussian_number,convolution_sigma):
    s = np.copy(t)
    mat = s
    if matrix_name != "raw":
    
        print "Normalizing with "+str(normalization)+" norm..."
        
        if normalization == "fragment-wise":
            floatorder = np.float64(order)
            s_norm_x = np.linalg.norm(s, ord=floatorder, axis=0)
            s_norm_y = np.linalg.norm(s, ord=floatorder, axis=1)
            s_norm = np.tensordot(s_norm_x,s_norm_y,axes=0)
            s[s_norm!=0] = s[s_norm!=0]/s_norm[s_norm!=0]
            print "Normalized "+str(normalization)+" with order "+str(order)
            
        elif normalization == "matrix-wise":
            
            floatorder = np.float64(order)
            s_norm = np.linalg.norm(s, ord=floatorder)
            s = s/s_norm
            print "Normalized "+str(normalization)+" with order "+str(order)
            
        elif normalization == "SCN":
            
            for iteration in range(1,iterations):
                sumrow = s.sum(axis=1)[:,None]
                sumcols = s.sum(axis=0)[None,:]
                s[sumrow!=0] = s[sumrow!=0]/sumrow[sumrow!=0]
                s[sumcols!=0] = s[sumcols!=0]/sumcols[sumcols!=0]
                print "Normalized "+str(iteration+1)+" time"+str("" if iteration <= 1 else "s")
            
            s = (s+s.T)/2
        
        elif normalization == "mirnylib":
            
            s_mirny = ntls.iterativeCorrection(s, iterations)[0]
            s = s_mirny
            print "Normalized "+str(iterations)+" time"+str("" if iterations <= 1 else "s")

        elif normalization == "sparsity":
            M = s.sum()
            sums = s.sum(axis=0)
            C = [[sums[i]*sums[j] for i in range(len(sums))] for j in range(len(sums))]/M
            s_coverage = s
            s_coverage[C!=0] /= C[C!=0]
            s = s_coverage
            
            print "Normalized for "+str(normalization)
            
        else:
            print "Error in normalization, using matrix-wise by default"
            s_norm = np.linalg.norm(s)
            s /= s_norm
        
        #Apply log or power
        try:
            s_exp = s**exposant
            s = s_exp
            print "Applied "+str(exposant)+" power to matrix"
        except ValueError:
            if exposant in ["log10", "log", "ln10"]:
                s = log10(s.astype(float))
                print "Applied base-10 logarithm to matrix"
            elif exposant in ["ln", "logarithm", "logarithme"]:
                s = log(s.astype(float))
                print "Applied natural logarithm to matrix"
            elif exposant in ["ln2", "log2"]:
                s = log2(s.astype(float))
                print "Applied base-2 logarithm to matrix"
            else:
                print "Warning, no valid normalization function encounter, ignoring"
        
        if matrix_name != "normalized":
            
            if "correlation" in matrix_name:
                s_corr = corrcoef(s)
                s_corr[s_corr<0] = 0
                s = s_corr
                print "Applied correlation function"
            
            if matrix_name != "correlation":
                
                if not "convolution" in matrix_name:
                    print "Error in matrix mode, using raw by default"
                    s = mat
                    
                else:
                    print "Convoluting..."
                    for i in range(0,gaussian_number):
                        s_gauss = ndimage.filters.gaussian_filter(s,convolution_sigma)
                        s = s_gauss
                        print "Convoluted "+str(i+1)+" time"+str("" if i+1 <= 1 else "s")
    return s

def square(s):
    n = np.min(s.shape)
    return s[0:n,0:n]

def manual_bin(s,factor=3):
    n,m = s.shape
    p = factor-1
    M = np.array([[np.sum(s[i:np.min([i+p,n]),j:np.min([j+p,m])]) for i in range(0,n-1,p+1)] for j in range(0,m-1,p+1)])
    return M
    
def directional(A, nw):
    
    n1 = A.shape[0]
    print "Size of the matrix entetered for the directional index:"
    print n1;
    signal1 = np.zeros((n1, 1));
    
    for i in range(0,n1) :
        vect_left = [];
        vect_right = [];
        
        for k in range(i-1,i-nw-1,-1) :
            kp =k;
            if k < 0 :
                kp = n1 +k ;
            if A[i,kp] > 0 :
                vect_left.append(math.log(A[i,kp]));
            else :
                vect_left.append(0);
    
    
        for k in range(i+1,i+nw+1) :
            kp =k;
            if k >= n1 :
                kp = k - n1;
            if A[i,kp] > 0 :
                vect_right.append(math.log(A[i,kp]));
            else :
                vect_right.append(0);

        if sum(vect_left) != 0 and sum(vect_right) != 0 :
            signal1[i] =  stats.ttest_rel(vect_right,vect_left)[0];
        else :
            signal1[i] =  0;
    
    return signal1
    
def trim_for_sparsity(M, n_std=3, s_min=None, s_max=None):
    sparsity = np.sum(M,0)
    mean = np.mean(sparsity)
    std = np.std(sparsity)
    if s_min is None:
        s_min = mean-n_std*std
    if s_max is None:
        s_max = mean+n_std*std
    
    f = (sparsity>s_min)*(sparsity<s_max)
    N = M[f][:,f]
    print N.shape
    return N

def rebin_kb(M,frags,kb=10):
    indices = []
    lengths = frags['len_bp']
    contigs = frags['id_c']
    
    for i in range(len(lengths)):
        if len(indices)==0 or np.sum(lengths[indices[-1]:i])>kb*1000.0 or contigs[indices[-1]]!=contigs[i]:
            indices.append(i)
    
    print indices
    row_bins = []           
    for j in range(len(indices)):
        print j
        print indices[j],indices[min(j+1,len(indices)-1)]
        rows_to_merge = M[indices[j]:indices[min(j+1,len(indices)-1)]]
        print rows_to_merge.shape
        row_bins.append(np.sum(rows_to_merge,axis=0))
    row_bins=np.array(row_bins)
    print row_bins
    print row_bins.shape
    N = []
    for j in range(len(indices)):
        columns_to_merge = row_bins[:,indices[j]:indices[min(j+1,len(indices)-1)]]
        print columns_to_merge
        print columns_to_merge.shape
        N.append(np.sum(columns_to_merge,axis=1))
    N = np.array(N)        
    print N.shape  
    return N
