import numpy as np
from numpy import arange, dot

# .\unwrapLaplacian.m

# Unwrap phase using Laplacian operation
# Schofield and Zhu, 2003 M.A. Schofield and Y. Zhu, Fast phase unwrapping
# algorithm for interferometric applications, Opt. Lett. 28 (2003), pp. 1194?196
#   [iFreq ] = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size)
# 
#   output
#   iFreq - unwrapped phase

#   input
#   iFreq_raw - A wrapped field map
#   matrix_size - dimension of the field of view
#   voxel_size - size of voxel in mm
#   
#   Created by Tian Liu on 2008.8.10
#   Modified by Tian Liu on 2011.01.26
#   Last modified by Tian Liu on 2013.07.24
#   Last modified by Tian Liu on 2014.02.10


# def del2(A,hx=1,hy=1,hz=1):
# 	L=A.copy()
# 	L=L.astype(np.float16)
# 	if len(A.shape)==2:
# 		for i in range (0,len(A)):
# 			for j in range (0,len(A[0])):
# 				if i==0:
# 					y=(2*A[i,j]+4*A[i+2,j]-5*A[i+1,j]-A[i+3,j])/4/hy/hy
# 				elif i==len(A)-1:
# 					y=(2*A[i,j]+4*A[i-2,j]-5*A[i-1,j]-A[i-3,j])/4/hy/hy
# 				else:
# 					y=(A[i+1,j]-2*A[i,j]+A[i-1,j])/4/hy/hy
# 				if j==0:
# 					x=(2*A[i,j]+4*A[i,j+2]-5*A[i,j+1]-A[i,j+3])/4/hx/hx
# 				elif j==len(A[0])-1:
# 					x=(2*A[i,j]+4*A[i,j-2]-5*A[i,j-1]-A[i,j-3])/4/hx/hx
# 				else:
# 					x=(A[i,j+1]-2*A[i,j]+A[i,j-1])/4/hx/hx
# 				L[i,j]=x+y
# 		return L
# 	elif len(A.shape)==3:
# 		for i in range (len(A)):
# 			for j in range (len(A[0])):
# 				for k in range (len(A[0,0])):
# 					if i==0:
# 						y=(2*A[i,j,k]+4*A[i+2,j,k]-5*A[i+1,j,k]-A[i+3,j,k])/6/hy/hy
# 					elif i==len(A)-1:
# 						y=(2*A[i,j,k]+4*A[i-2,j,k]-5*A[i-1,j,k]-A[i-3,j,k])/6/hy/hy
# 					else:
# 						y=(A[i+1,j,k]-2*A[i,j,k]+A[i-1,j,k])/6/hy/hy
# 					if j==0:
# 						x=(2*A[i,j,k]+4*A[i,j+2,k]-5*A[i,j+1,k]-A[i,j+3,k])/6/hx/hx
# 					elif j==len(A[0])-1:
# 						x=(2*A[i,j,k]+4*A[i,j-2,k]-5*A[i,j-1,k]-A[i,j-3,k])/6/hx/hx
# 					else:
# 						x=(A[i,j+1,k]-2*A[i,j,k]+A[i,j-1,k])/6/hz/hz
# 					if k==0:
# 						z=(2*A[i,j,k]+4*A[i,j,k+2]-5*A[i,j,k+1]-A[i,j,k+3])/6/hz/hz
# 					elif k==len(A[0,0])-1:
# 						z=(2*A[i,j,k]+4*A[i,j,k-2]-5*A[i,j,k-1]-A[i,j,k-3])/6/hz/hz
# 					else:
# 						z=(A[i,j,k+1]-2*A[i,j,k]+A[i,j,k-1])/6/hz/hz
# 					L[i,j,k]=x+y+z
# 		return L


def del2(f, voxel_size):
# %DEL2 Discrete Laplacian.
# %   L = DEL2(U), when U is a matrix, is a discrete approximation of
# %   0.25*del^2 u = (d^2u/dx^2 + d^2u/dy^2)/4.  The matrix L is the same
# %   size as U, with each element equal to the difference between an 
# %   element of U and the average of its four neighbors.
# %
# %   L = DEL2(U), when U is an N-D array, returns an approximation of
# %   (del^2 u)/2/n, where n is ndims(u).
# %
# %   L = DEL2(U,H), where H is a scalar, uses H as the spacing between
# %   points in each direction (H=1 by default).
# %
# %   L = DEL2(U,HX,HY), when U is 2-D, uses the spacing specified by HX
# %   and HY. If HX is a scalar, it gives the spacing between points in
# %   the x-direction. If HX is a vector, it must be of length SIZE(U,2)
# %   and specifies the x-coordinates of the points.  Similarly, if HY
# %   is a scalar, it gives the spacing between points in the
# %   y-direction. If HY is a vector, it must be of length SIZE(U,1) and
# %   specifies the y-coordinates of the points.
# %
# %   L = DEL2(U,HX,HY,HZ,...), when U is N-D, uses the spacing given by
# %   HX, HY, HZ, etc. 
# %
# %   Class support for input U:
# %      float: double, single
# %
# %   See also GRADIENT, DIFF.
# 
# %   Copyright 1984-2015 The MathWorks, Inc. 
    v = voxel_size.copy()
    f, ndim, loc, rflag = parse_input(f, v)
    
    # Loop over each dimension. Permute so that the del2 is always taken along
    # the columns.
    if ndim == 1:
        perm = [0, 1];
    else:
        perm = np.arange(1,ndim+1) # Cyclic permutation
        perm[-1]=0
        
    v = np.zeros(f.shape, f.dtype)
    for k in range(ndim):
        n = f.shape[0]
        x = loc[k].flatten()
        h = np.diff(x)
        g = np.zeros(f.shape, f.dtype) # % case of singleton dimension

        # Take centered second differences on interior points to compute g/2
        if n > 2:
            h_repeat = np.repeat(h[:,np.newaxis],f.shape[1],axis=1)
            if n>3: h_repeat = np.repeat(h_repeat[...,np.newaxis],f.shape[2],axis=2)
            g[1:n-1,:] = (np.diff(f[1:n,:],axis=0) / h_repeat[1:n-1,:] - np.diff(f[0:n-1,:],axis=0) / h_repeat[0:n-2,:]) / (h_repeat[1:n-1,:] + h_repeat[0:n-2,:])
                                                            
        # Linearly extrapolate second differences from interior
        if n > 3:
            g[0,:]   =  g[1,:] * (h[0]+h[1])/h[1] - g[2,:]*(h[1])/h[2]
            g[n-1,:] = -g[n-3,:]*(h[n-2]) / h[n-3] + g[n-2,:] * (h[n-2]+h[n-3])/h[n-3]
        elif n==3:
            g[0,:] = g[1,:]
            g[n-1,:] = g[1,:]
        else:
            g[0,:]=0
            g[n-1,:]=0
    
        if ndim==1:
            v = v + g
        else:
            order = np.concatenate((np.arange(k,ndim), np.arange(0,k)),axis=0)
            iorder = np.zeros(order.shape)
            iorder[order] = np.arange(3).astype(int)
            v = v + np.transpose(g,iorder.astype(int))
    
        # Set up for next pass through the loop
        f = np.transpose(f,perm)
    
    v = v / ndim
    if rflag: 
        v = np.transpose(v) 
    return v   
def parse_input(f, v):
    
    # Flag vector case and column vector case.
    ndim = len(f.shape);
    loc = [[]*ndim]
    vflag = 0
    rflag = 0
    
    indx = f.shape
    
    if len(f.shape) == len(v): # del2(f,hx,hy,hz,...)
        # Swap 1 and 2 since x is the second dimension and y is the first.
        loc = v.tolist()
        if ndim >1:
            loc[:2] = loc[:2][::-1]
        
        # replace any scalar step-size with corresponding position vector
        for k in range(ndim):
            loc[k] = loc[k]*np.arange(1,indx[k]+1)
    return f, ndim, loc, rflag

def unwrapLaplacian(iFreq_raw=None,matrix_size=None,voxel_size=[1,1,1]):
    if len(matrix_size) == 2:
        matrix_size[2]=1

    if len(voxel_size) == 2:
        voxel_size[2]=1

    Y,X,Z=np.meshgrid(arange(- matrix_size[1] / 2,matrix_size[1] / 2 ),arange(- matrix_size[0] / 2,matrix_size[0] / 2 ),arange(- matrix_size[2] / 2,matrix_size[2] / 2 ))
    X=dot(X,voxel_size[0])
    Y=dot(Y,voxel_size[1])
    Z=dot(Z,voxel_size[2])

    if matrix_size[2] > 1:
        h=((X == 0) * (Y == 0) * (Z == 0)).astype(np.float32) #np.array((logical_and(logical_and((X == 0),(Y == 0)),(Z == 0)))*1)
        k=6 * del2(h,voxel_size)
        kernel=np.fft.fftn(np.fft.fftshift(k))
    else:
        h=((X == 0) * (Y == 0)).astype(np.float16)#(logical_and((X == 0),(Y == 0)))*1
        k=4 * del2(h,voxel_size[0:2])
        kernel=np.fft.fftn(np.fft.fftshift(k))

    inv_kernel=1.0 / kernel
    inv_kernel[np.isinf(inv_kernel)]=0
    inv_kernel[abs(kernel) < 1e-10]=0
    first_term=np.cos(iFreq_raw) * np.fft.ifftn(kernel*np.fft.fftn(np.sin(iFreq_raw)))

    second_term=np.sin(iFreq_raw) * np.fft.ifftn(kernel * np.fft.fftn(np.cos(iFreq_raw)))

    phi_est=np.fft.ifftn(inv_kernel.astype(float) * np.fft.fftn((first_term.astype(np.float32) - second_term.astype(np.float32))))
    return phi_est.astype(np.float32)

# if __name__=='__main__':

# 	A=np.array([[3,4,6,7],[8,9,100,11],[12,13,14,15],[16,17,18,19]])
# 	print(A)
# 	print(del2(A))
# 	print(del2(A,2,1))
# 	A=np.array([[[0 for i in range(3)]for i in range(1)]for i in range (4)])
# 	A[:,:,0]=np.array([[3],[8],[12],[16]])
# 	A[:,:,1]=A[:,:,0].copy()
# 	A[:,:,2]=A[:,:,1].copy()
# 	print(len(A))
# 	print(len(A[0]))
# 	print(len(A[0,0]))
# 	print('A=\n',A,'\n')
# 	print(A[1,2,3])
# 	print(del2(A)[:,:,0])