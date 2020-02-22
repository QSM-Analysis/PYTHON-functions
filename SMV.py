# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
from numpy import dot,max
from sphere_kernel import*
# .\SMV.m

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
    

def SMV(iFreq=None,varargin=None,*args,**kwargs):
    if 1 == len(varargin):
        K=varargin[0]
# .\SMV.m:19
    else:
        matrix_size=varargin[0]
# .\SMV.m:21
        voxel_size=varargin[1]
# .\SMV.m:22
        if (len(varargin) < 3):
            radius=dot(round(6 / max(voxel_size)),max(voxel_size))
# .\SMV.m:24
        else:
            radius=varargin[2]
# .\SMV.m:26
        K=sphere_kernel(matrix_size,voxel_size,radius)
# .\SMV.m:28
    
    y=np.fft.ifftn(multiply(np.fft.fftn(iFreq),K))
# .\SMV.m:31
    return y