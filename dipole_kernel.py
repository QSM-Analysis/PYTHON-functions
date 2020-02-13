# Generated with SMOP  0.41
import numpy as np
from smop.libsmop import *
import os
import sys
import scipy.io as scio
def dipole_kernel(varargin=None,*args,**kwargs):
    
    matrix_size,voxel_size,Bo_dir,domain = parse_inputs(varargin[:])
    
    matrix_size=scio.loadmat('matlabdata_matrix_size_FOR RDF.mat')['matrix_size']
    voxel_size=scio.loadmat('matlabdata_voxel_size_FOR RDF.mat')['voxel_size']
    B0_dir=scio.loadmat('matlabdata_ B0_dir_FOR RDF.mat')['B0_dir']
    matrix_size=ravel(matrix_size)
    
    if Bo_dir == 1:
        item = np.array([1,0,0])
        Bo_dir = item.conj().T
    elif Bo_dir == 2:
        Bo_dir = np.array([0,1,0]).conj().T
    elif Bo_dir == 3:
        Bo_dir = np.array([0,0,1]).conj().T
    print('B0_dir=')
    print(B0_dir)
    temp1=np.zeros((64,1))
    temp1[0]=B0_dir[0]
    temp2=np.zeros((64,1))        
    temp2[1]=B0_dir[1]
    temp3=np.zeros((64,1))
    temp3[2]=B0_dir[2]
    if domain=='kspace':
        (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,matrix_size[1]/2), np.arange(-matrix_size[0]/2,matrix_size[0]/2), np.arange(-matrix_size[2]/2,matrix_size[2]/2))
# .\dipole_kernel.m:39
        X=X / (matrix_size[0]*voxel_size[0])
# .\dipole_kernel.m:43
        Y=Y / (matrix_size[1]*voxel_size[1])
# .\dipole_kernel.m:44
        Z=Z / (matrix_size[2]*voxel_size[2])
# .\dipole_kernel.m:45
        D=1 / 3 - (dot(X,temp1) + dot(Y,temp2) + dot(Z,temp3)) ** 2.0 / (X ** 2 + Y ** 2 + Z ** 2)
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
            Y=dot(Y,voxel_size[1])
# .\dipole_kernel.m:57
            Z=dot(Z,voxel_size[2])
            
# .\dipole_kernel.m:58
            d=(dot(3,(dot(X,temp1) + dot(Y,temp2) + dot(Z,temp3)) ** 2) - X ** 2 - Y ** 2 - Z ** 2) / (dot(dot(4,pi),(X ** 2 + Y ** 2 + Z ** 2) ** 2.5))
# .\dipole_kernel.m:60
            d[np.isnan(d)]=0
# .\dipole_kernel.m:62
            D=np.fft.fftn(np.fft.fftshift(d))
# .\dipole_kernel.m:63
    
    return D

#def parse_inputs(varargin=None,*args,**kwargs):
#    if np.size(varargin) < 3:
##        error('At least matrix_size, voxel_size and B0_dir are required')
#        tkinter.messagebox.showwarning('警告','At least matrix_size, voxel_size and B0_dir are required')        
#    print('varargin=',varargin)
#    print('args=',args)
#    print('kwargs=',kwargs)
#    matrix_size=varargin[0]
## .\dipole_kernel.m:76
#    voxel_size=varargin[1]
## .\dipole_kernel.m:77
#    B0_dir=varargin[2]
## .\dipole_kernel.m:78
#    domain='kspace'
## .\dipole_kernel.m:79
#    if np.size(varargin) > 3:
#        for k in range(3,np.size(varargin)):
#            if varargin[k].lower()=='imagespace':
#                domain='imagespace'
## .\dipole_kernel.m:84 
#    return matrix_size,voxel_size,B0_dir,domain
    
def parse_inputs(varargin):
    if np.size(varargin) < 3:
        print("At least matrix_size, voxel_size and B0_dir are required")
        os._exit(0)
    print('varargin=',varargin)
    [[256,256,64],]
    matrix_size = varargin[0]
    print(matrix_size)
    voxel_size = varargin[1]
    Bo_dir = varargin[2]
    domain = 'kspace'

    if varargin.shape[0] > 3:
        for k in range(3,varargin.shape[0]):
            if varargin[k] == 'imagespace':
                domain = 'imagespace'
    return [matrix_size,voxel_size,Bo_dir,domain]