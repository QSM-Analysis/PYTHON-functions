# -*- coding: utf-8 -*-
import numpy as np
from numpy import dot,max, multiply
from PyQSM.sphere_kernel import sphere_kernel

# Spherical Mean Value operator
#   y=SMV(iFreq,matrix_size,voxel_size,radius)
#   
#   output
#   y - reulstant image after SMV
# 
#   input
#   iFreq - input image
#   matrix_size - dimension of the field of view
#   voxel_size - the size of the voxel
#   radius - radius of the sphere in mm
    
    #   Created by Tian Liu in 2010
#   Last modified by Tian Liu on 2013.07.24
    

def SMV(iFreq=None,varargin=None):

    if 1 == len(varargin):
        K=varargin[0]
    else:
        matrix_size=varargin[0]
        voxel_size=varargin[1]
        if (len(varargin) < 3):
            radius=dot(round(6 / max(voxel_size)),max(voxel_size))
        else:
            radius=varargin[2]
        K = sphere_kernel(matrix_size,voxel_size,radius)
    
    y=np.fft.ifftn(np.fft.fftn(iFreq) *K)
    return y