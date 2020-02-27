# -*- coding: utf-8 -*-
# Created on Feb 27, 2020
# Author: Maxwell4444
# All right reserved

# r2 = arlo(te,y)
#    
# output
#   r2 r2-map (in Hz)
#      empty when only one echo is provided
# 
# input
#   te  array containing te values (in s)
#   y   a multi-echo data set of arbitrary dimension
#       echo should be the last dimension
#       
# If you use this function please cite
# 
# Pei M, Nguyen TD, Thimmappa ND, Salustri C, Dong F, Cooper MA, Li J, 
# Prince MR, Wang Y. Algorithm for fast monoexponential fitting based 
# on Auto-Regression on Linear Operations (ARLO) of data. 
# Magn Reson Med. 2015 Feb;73(2):843-50. doi: 10.1002/mrm.25137. 
# Epub 2014 Mar 24. PubMed PMID: 24664497; 
# PubMed Central PMCID:PMC4175304.

import numpy as np
_eps = 0.000000001

def arlo(te, y):
    nte = te.size
    if nte < 2: return []
    
    sz = y.shape
    edx = len(sz)
    if sz[edx-1]!= nte:
        print('Error for nte')
        
    yy  = np.zeros(sz[:-1])
    yx  = np.zeros(sz[:-1])
    beta_yx = np.zeros(sz[:-1])
    beta_xx = np.zeros(sz[:-1])
    sl = []
    dl = []
    
    for j in range(nte-2):
        alpha = (te[j+2]-te[j])*(te[j+2]-te[j])/2/(te[j+1]-te[j])
        tmp = (2*te[j+2]*te[j+2] - te[j]*te[j+2] - te[j]*te[j] + 3*te[j]*te[j+1] -3*te[j+1]*te[j+2])/6 
        beta = tmp/(te[j+2]-te[j+1])
        gamma = tmp/(te[j+1]-te[j])
        
        y1 = y[...,j]*(te[j+2]-te[j]-alpha+gamma)+y[...,j+1]*(alpha-beta-gamma)+y[...,j+2]*beta
        x1 = y[...,j]-y[...,j+2]
        yy = yy + y1 * y1
        yx = yx + y1 * x1
        beta_yx = beta_yx + beta*y1*x1
        beta_xx = beta_xx + beta*x1*x1
    
    r2 = (yx + beta_xx)/(beta_yx + yy + _eps)
    r2[np.isnan(r2)] = 0
    r2[np.isinf(r2)] = 0
    
    return r2