# Generated with SMOP  0.41
import numpy as np

import os
import sys
import scipy.io as scio
def dipole_kernel(matrix_size,voxel_size,B0_dir,domain='kspace'):
    #matrix_size,voxel_size,B0_dir,domain = parse_inputs(varargin[:])
    #matrix_size=matrix_size.flatten(1)
    if B0_dir[0] == 1:
        B0_dir = np.array([1,0,0]).conj().T
        
    elif B0_dir[1] == 1:
        B0_dir = np.array([0,1,0]).conj().T
    elif B0_dir[2] == 1:
        B0_dir = np.array([0,0,1]).conj().T
    if domain=='kspace':
        (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,matrix_size[1]/2), np.arange(-matrix_size[0]/2,matrix_size[0]/2), np.arange(-matrix_size[2]/2,matrix_size[2]/2))

        X=X / (matrix_size[0]*voxel_size[0])
        Y=Y / (matrix_size[1]*voxel_size[1])
        Z=Z / (matrix_size[2]*voxel_size[2])

        D=1 / 3 - (X*B0_dir[0] + Y*B0_dir[1] +Z*B0_dir[2]) ** 2.0 / (X ** 2 + Y ** 2 + Z ** 2)
        D[np.isnan(D)]=0
        D=np.fft.fftshift(D)
    else:
        if domain=='imagespace':
            (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,(matrix_size[1]/2-1)), 
                                np.arange(-matrix_size[0]/2,(matrix_size[0]/2-1)), 
                                np.arange(-matrix_size[2]/2,(matrix_size[2]/2-1)) )
            
            X=X*voxel_size[0]
            Y=Y*voxel_size[1]
            Z=Z*voxel_size[2]
            
            d=(3*(X*B0_dir[0] + Y*B0_dir[1] +Z*B0_dir[2]) ** 2 - X ** 2 - Y ** 2 - Z ** 2) / (4 * np.pi * (X ** 2 + Y ** 2 + Z ** 2) ** 2.5)
            
            d[np.isnan(d)]=0
            D=np.fft.fftn(np.fft.fftshift(d)).astype(float)
    
    return D
    
# def parse_inputs(varargin):
#     if len(varargin) < 3:
#         print("At least matrix_size, voxel_size and B0_dir are required")
#         os._exit(0)
#     matrix_size = varargin[0]
#     voxel_size = varargin[1]
#     B0_dir = varargin[2]
#     
#     domain = 'kspace'
#     if len(varargin) > 3:
#         for k in range(3,len(varargin)):
#             if varargin[k] == 'imagespace':
#                 domain = 'imagespace'
#     return [matrix_size,voxel_size,B0_dir,domain]