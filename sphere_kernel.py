# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
from numpy import dot,sum
# .\sphere_kernel.m

    # Generate a Spherical kernel with the sum normalized to one
#   y = SMV_kernel(matrix_size,voxel_size, radius)
#   
#   output
#   y - kernel
# 
#   input
#   matrix_size - the dimension of the field of view
#   voxel_size - the size of the voxel in mm
#   radius - the raidus of the sphere in mm
    
    #   Created by Tian Liu in 2010
#   Modified by Tian on 2011.02.01
#   Modified by Tian on 2011.03.14 The sphere is now rendered.
#   Last modified by Tian Liu on 2013.07.23
  
def sphere_kernel(matrix_size=None,voxel_size=None,radius=None,*args,**kwargs):
    radius=float(radius)
    matrix_size=np.array(matrix_size[0])
# .\dipole_kernel.m:76
    #voxel_size=varargin[1]
    voxel_size=np.array(voxel_size.T[0])
    (Y,X,Z)=np.meshgrid(np.arange(-matrix_size[1]/2,matrix_size[1]/2), np.arange(-matrix_size[0]/2,matrix_size[0]/2), np.arange(-matrix_size[2]/2,matrix_size[2]/2))
# .\sphere_kernel.m:19
    X=dot(X,voxel_size[0])
# .\sphere_kernel.m:23
    Y=dot(Y,voxel_size[1])
# .\sphere_kernel.m:24
    Z=dot(Z,voxel_size[2])
# .\sphere_kernel.m:25
    a=abs(X) - dot(0.5,voxel_size[0])
    a[a<0]=0
    b=abs(Y) - dot(0.5,voxel_size[1])
    b[b<0]=0
    c=abs(Z) - dot(0.5,voxel_size[2])
    c[c<0]=0
    Sphere_out=(a ** 2 + b ** 2 + c ** 2) > (radius ** 2)
# .\sphere_kernel.m:26
    Sphere_in=((abs(X) + dot(0.5,voxel_size[0])) ** 2 + (abs(Y) + dot(0.5,voxel_size[1])) ** 2 + (abs(Z) + dot(0.5,voxel_size[2])) ** 2) <= radius ** 2
# .\sphere_kernel.m:30
    Sphere_mid=np.zeros(matrix_size.astype(int))
# .\sphere_kernel.m:35
    split=10
# .\sphere_kernel.m:37
    
    (X_v,Y_v,Z_v)=np.meshgrid(np.arange(- split + 0.5,split + 0.5),np.arange(- split + 0.5,split + 0.5),np.arange(- split + 0.5,split + 0.5))
# .\sphere_kernel.m:38
    X_v=X_v / (dot(2,split))
# .\sphere_kernel.m:39
    Y_v=Y_v / (dot(2,split))
# .\sphere_kernel.m:40
    Z_v=Z_v / (dot(2,split))
# .\sphere_kernel.m:41
    shell=1 - Sphere_in - Sphere_out
# .\sphere_kernel.m:43
    X=X[shell == 1]
# .\sphere_kernel.m:44
    Y=Y[shell == 1]
# .\sphere_kernel.m:45
    Z=Z[shell == 1]
# .\sphere_kernel.m:46
    shell_val=np.zeros(X.shape)
# .\sphere_kernel.m:47
    for i in range(len(X)):
        xx=X[i]
# .\sphere_kernel.m:50
        yy=Y[i]
# .\sphere_kernel.m:51
        zz=Z[i]
# .\sphere_kernel.m:52
        occupied=((xx + dot(X_v,voxel_size[0])) ** 2 + (yy + dot(Y_v,voxel_size[1])) ** 2 + (zz + dot(Z_v,voxel_size[2])) ** 2) <= radius ** 2
# .\sphere_kernel.m:54
        shell_val[i]=sum(occupied[:]) / np.size(X_v)
# .\sphere_kernel.m:57
    
    Sphere_mid[shell == 1]=shell_val
# .\sphere_kernel.m:60
    Sphere=Sphere_in + Sphere_mid
# .\sphere_kernel.m:62
    Sphere=Sphere / sum(Sphere[:])
# .\sphere_kernel.m:63
    y=np.fft.fftn(np.fft.fftshift(Sphere))
# .\sphere_kernel.m:64
    return y

if __name__=='__main__':
    pass