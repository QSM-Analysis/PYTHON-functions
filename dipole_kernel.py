# Generated with SMOP  0.41
<<<<<<< HEAD
import numpy as np
#from smop.libsmop import *
import os
import sys
import scipy.io as scio
from numpy import ravel, dot
=======
from smop.libsmop import *
import numpy as np
import tkinter.messagebox
from numpy import dot,multiply,pi
# .\dipole_kernel.m

    # Generation of the Dipole Kernel
#   D = dipole_kernel(matrix_size, voxel_size, B0_dir, ...)
    
    #   output
#   D - dipole kernel saved in Fourier space
# 
#   input
#   matrix_size - the size of the matrix
#   voxel_size - the size of a voxel
#   B0_dir - the direction of the B0 field
#   domain - 'imagespace' or 'kspace'
#       Fourier domain expression:
#       Salomir et al. Concepts in Magn Reson Part B 2003, 19(B):26-34
#       Image domain expression: 
#       Li et al. Magn Reson Med 2004, 51(5):1077-82
    
    
    #   Created by Tian Liu in 2008
#   Modified by Tian Liu on 2011.02.01
#   Last modified by Tian Liu on 2013.07.22
#   
def parse_inputs(varargin=None,*args,**kwargs):
    if varargin.shape[0] < 3:
#        error('At least matrix_size, voxel_size and B0_dir are required')
        tkinter.messagebox.showwarning('警告','At least matrix_size, voxel_size and B0_dir are required') 
    
    matrix_size=np.array(varargin[0][0])
# .\dipole_kernel.m:76
    #voxel_size=varargin[1]
    voxel_size=np.array((varargin[1]).T[0])
# .\dipole_kernel.m:77
    B0_dir=varargin[2]
# .\dipole_kernel.m:78
    domain='kspace'
# .\dipole_kernel.m:79
    if varargin.shape[0] > 3:
        for k in range(3,varargin.shape[0]):
            if varargin[k].lower()=='imagespace':
                domain='imagespace'
# .\dipole_kernel.m:84 
    return matrix_size,voxel_size,B0_dir,domain
    
>>>>>>> a4332f83f492ca65591f7d1ad49b58f48c7db3bb
def dipole_kernel(varargin=None,*args,**kwargs):
    matrix_size,voxel_size,B0_dir,domain=parse_inputs(varargin[:])
# .\dipole_kernel.m:28
    if (B0_dir == '1'):
        B0_dir=np.array([1,0,0])
# .\dipole_kernel.m:31
    else:
        if (B0_dir == '2'):
            B0_dir=np.array([0,1,0])
# .\dipole_kernel.m:33
        else:
            if (B0_dir == '3'):
                B0_dir=np.array([0,0,1])
# .\dipole_kernel.m:35
    if domain=='kspace':
        (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,matrix_size[1]/2), np.arange(-matrix_size[0]/2,matrix_size[0]/2), np.arange(-matrix_size[2]/2,matrix_size[2]/2))
# .\dipole_kernel.m:39
        X=X / (matrix_size[0]*voxel_size[0])
# .\dipole_kernel.m:43
        Y=Y / (matrix_size[1]*voxel_size[1])
# .\dipole_kernel.m:44
        Z=Z / (matrix_size[2]*voxel_size[2])

# .\dipole_kernel.m:45
        D=1 / 3 - (dot(X,B0_dir[0]) + dot(Y,B0_dir[1]) + dot(Z,B0_dir[2])) ** 2.0 / (X ** 2 + Y ** 2 + Z ** 2)
# .\dipole_kernel.m:47
        D[np.isnan(D)]=0
# .\dipole_kernel.m:48
        D=np.fft.fftshift(D)
# .\dipole_kernel.m:49
    else:
        if domain=='imagespace':
            (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,matrix_size[1]/2), np.arange(-matrix_size[0]/2,matrix_size[0]/2), np.arange(-matrix_size[2]/2,matrix_size[2]/2))
# .\dipole_kernel.m:52
            X=dot(X,voxel_size[0])
# .\dipole_kernel.m:56
            Y=dot(Y,voxel_size[1])
# .\dipole_kernel.m:57
            Z=dot(Z,voxel_size[2])
# .\dipole_kernel.m:58
            d=(dot(3,(dot(X,B0_dir[0]) + dot(Y,B0_dir[1]) + dot(Z,B0_dir[2])) ** 2) - X ** 2 - Y ** 2 - Z ** 2) / (dot(dot(4,pi),(X ** 2 + Y ** 2 + Z ** 2) ** 2.5))
# .\dipole_kernel.m:60
            d[np.isnan(d)]=0
# .\dipole_kernel.m:62
            D=np.fft.fftn(np.fft.fftshift(d))
# .\dipole_kernel.m:63
    
    return D
    
if __name__ == '__main__':
    pass
    