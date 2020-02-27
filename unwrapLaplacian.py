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

class opts0():
    def __init__(self):
        self.matrix_size=None
        self.voxel_size = np.array([1,1,1])
defopts=opts0()

def del2(A,hx=1,hy=1,hz=1):
	L=A.copy()
	L=L.astype(np.float16)
	if len(A.shape)==2:
		for i in range (0,len(A)):
			for j in range (0,len(A[0])):
				if i==0:
					y=(2*A[i,j]+4*A[i+2,j]-5*A[i+1,j]-A[i+3,j])/4/hy/hy
				elif i==len(A)-1:
					y=(2*A[i,j]+4*A[i-2,j]-5*A[i-1,j]-A[i-3,j])/4/hy/hy
				else:
					y=(A[i+1,j]-2*A[i,j]+A[i-1,j])/4/hy/hy
				if j==0:
					x=(2*A[i,j]+4*A[i,j+2]-5*A[i,j+1]-A[i,j+3])/4/hx/hx
				elif j==len(A[0])-1:
					x=(2*A[i,j]+4*A[i,j-2]-5*A[i,j-1]-A[i,j-3])/4/hx/hx
				else:
					x=(A[i,j+1]-2*A[i,j]+A[i,j-1])/4/hx/hx
				L[i,j]=x+y
		return L
	elif len(A.shape)==3:
		for i in range (len(A)):
			for j in range (len(A[0])):
				for k in range (len(A[0,0])):
					if i==0:
						y=(2*A[i,j,k]+4*A[i+2,j,k]-5*A[i+1,j,k]-A[i+3,j,k])/6/hy/hy
					elif i==len(A)-1:
						y=(2*A[i,j,k]+4*A[i-2,j,k]-5*A[i-1,j,k]-A[i-3,j,k])/6/hy/hy
					else:
						y=(A[i+1,j,k]-2*A[i,j,k]+A[i-1,j,k])/6/hy/hy
					if j==0:
						x=(2*A[i,j,k]+4*A[i,j+2,k]-5*A[i,j+1,k]-A[i,j+3,k])/6/hx/hx
					elif j==len(A[0])-1:
						x=(2*A[i,j,k]+4*A[i,j-2,k]-5*A[i,j-1,k]-A[i,j-3,k])/6/hx/hx
					else:
						x=(A[i,j+1,k]-2*A[i,j,k]+A[i,j-1,k])/6/hz/hz
					if k==0:
						z=(2*A[i,j,k]+4*A[i,j,k+2]-5*A[i,j,k+1]-A[i,j,k+3])/6/hz/hz
					elif k==len(A[0,0])-1:
						z=(2*A[i,j,k]+4*A[i,j,k-2]-5*A[i,j,k-1]-A[i,j,k-3])/6/hz/hz
					else:
						z=(A[i,j,k+1]-2*A[i,j,k]+A[i,j,k-1])/6/hz/hz
					L[i,j,k]=x+y+z
		return L


def unwrapLaplacian(iFreq_raw=None,matrix_size=None,voxel_size=defopts.voxel_size):
    if len(matrix_size) == 2:
        matrix_size[2]=1

    if len(voxel_size) == 2:
        voxel_size[2]=1

    Y,X,Z=np.meshgrid(arange(- matrix_size[1] / 2,matrix_size[1] / 2 ),arange(- matrix_size[0] / 2,matrix_size[0] / 2 ),arange(- matrix_size[2] / 2,matrix_size[2] / 2 ))
    X=dot(X,voxel_size[0])
    Y=dot(Y,voxel_size[1])
    Z=dot(Z,voxel_size[2])

    if matrix_size[2] > 1:
        h=((X == 0) * (Y == 0) * (Z == 0)).astype(np.float16) #np.array((logical_and(logical_and((X == 0),(Y == 0)),(Z == 0)))*1)
        k=6 * del2(h,voxel_size[0],voxel_size[1],voxel_size[2])
        kernel=np.fft.fftn(np.fft.fftshift(k))
    else:
        h=((X == 0) * (Y == 0)).astype(np.float16)#(logical_and((X == 0),(Y == 0)))*1
        k=4 * del2(h,voxel_size[0],voxel_size[1])
        kernel=np.fft.fftn(np.fft.fftshift(k))

    inv_kernel=1.0 / kernel
    inv_kernel[np.isinf(inv_kernel)]=0
    inv_kernel[abs(kernel) < 1e-10]=0
    first_term=np.cos(iFreq_raw) * np.fft.ifftn(kernel*np.fft.fftn(np.sin(iFreq_raw)))

    second_term=np.sin(iFreq_raw) * np.fft.ifftn(kernel * np.fft.fftn(np.cos(iFreq_raw)))

    phi_est=np.fft.ifftn(inv_kernel * np.fft.fftn((first_term - second_term)))
    return phi_est

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