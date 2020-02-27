# Generated with SMOP  0.41
# Fit_ppm_complex.m

    # Projection onto Dipole Fields (PDF)
#   [p1, dp1, relres, p0]=Fit_ppm_complex(M)
#    
#   output
#   p1 - field map, may need further unwrapping
#   dp1 - a priori error estimate
#   relres - relative residual
#   p0 - initial phase
    
    #   input
#   M - a multi-echo and could be a multi-channel dataset
#       echo needs to be the 4th dimension
#       channel needs to be the 5th dimension
    
    #   When using the code, please cite 
#   T. Liu et al. MRM 2013;69(2):467-76
#   B. Kressler et al. IEEE TMI 2010;29(2):273-81
#   de Rochefort et al. MRM 2008;60(4):1003-1009
    
    #   The coil combination method is similar to
#   MA. Bernstein et al. MRM 1994;32:330-334
    
    #   Adapted from a linear fitting created by Ludovic de Rochefort
#   Modified by Tian Liu on 2011.06.01
#   Last modified by Alexey Dimov on 2016.05.12

import sys    
import numpy as np
from numpy import zeros, min, pi, abs,ones, dot, multiply, exp, mod, conj, angle, sqrt,\
    arange, floor, ceil, sum
from numpy.matlib import repmat
# dytpe='single'
_esp = 0.00000001

class opts0:
    def __init__(self):
        self.reltol= 0.0001
        self.max_iter = 30
defopts=opts0()

def Fit_ppm_complex(M=None,opts=defopts):
#     defopts.reltol = copy(0.0001)
#     defopts.max_iter = copy(30)
    
    #bytes=M.dtype
    bytes=len(M)*65536*16
    
    if (bytes > 1000000000.0):
        sz = M.shape
        p1=zeros(M.shape[:3],M.dtype)
        dp1=zeros(M.shape[:3],M.dtype)
        relres=zeros(M.shape[:3],M.dtype)
        p0=zeros(M.shape[:3],M.dtype)
        n=ceil(bytes / 1000000000.0)
        ns=floor(sz(3) / n)
        for slc in arange(1,sz(3),ns).reshape(-1):
            rng=arange(slc,min(slc + ns - 1,sz(3)))
            print('fitting slice '+str(rng(1))+' through '+str(rng(-1),' ...'))
            p1[:,:,rng],dp1[:,:,rng],relres[:,:,rng],p0[:,:,rng]=Fit_ppm_complex(M[:,:,rng,:,:],nargout=4)
        
            
        return p1,dp1,relres,p0
    
    
    #Modification to handle one echo datasets - assuming zero phase at TE = 0;
    #- AD, 05/12/2016
    if M.shape[3] == 1:
    #   M = cat(4,abs(M),M)
        M = np.concatenate((abs(M),M),axis=3)

    if len(M.shape)>4:
        
        if M.shape[4] > 1:
        # combine multiple coils together, assuming the coil is the fifth dimension
            M=np.sum(multiply(M,conj(repmat(M[:,:,:,1,:],concat([1,1,1,shape(M,4),1])))),5)
            M=multiply(sqrt(abs(M)),exp(dot(1j,angle(M))))
    
    # determine angle
    M=np.conjugate(M)
    s0=M.shape
    L_s0=len(s0)
    nechos = M.shape[-1]
    
    M=M.reshape([np.prod(s0[:(L_s0-1)]),s0[L_s0-1]])
    s=M.shape
    
    Y=np.angle(M[:,:np.min([3,nechos])])
    c=Y[:,1] - Y[:,0]
    ind = np.argmin([abs(c-2*pi),abs(c),abs(c+2*pi)],axis=0)
    c[ind==0]=c[ind==0]-2*np.pi
    c[ind==2]=c[ind==2]+2*np.pi
    
    for n in range(min([2,nechos-1])):
        cd=((Y[:,n + 1] - Y[:,n])) - c
        Y[cd < - pi,(n + 1):] = Y[cd < - pi,(n + 1):] + 2*pi
        Y[cd >   pi,(n + 1):] = Y[cd >   pi,(n + 1):] + 2*pi
    A=np.array([[1,0],[1,1],[1,2]])
    ip=dot(np.mat(A[:min([3,nechos]),:]).I,Y[:,:min([3,nechos])].conj().T)
    ip=np.asarray(ip)
    
    p0=np.array(ip[0,:].conj().T)
    p1=np.array(ip[1,:].conj().T)
    p0, p1 = p0.reshape(p0.size,1), p1.reshape(p1.size,1)
    
    dp1=p1.copy()
    tol = dot(np.linalg.norm(p1[:]),opts.reltol)
    miter=0
    
    # weigthed least square
    # calculation of WA'*WA
    v1=ones((1,nechos),M.dtype)
    v2=np.arange(0,nechos).astype(M.dtype).reshape(1,nechos)
    a11=np.sum(multiply(abs(M) ** 2.0,dot(ones((s[0],1),M.dtype),v1 ** 2)),1)
    a12=np.sum(multiply(abs(M) ** 2.0,dot(ones((s[0],1),M.dtype),v1 * v2)),1)
    a22=np.sum(multiply(abs(M) ** 2.0,dot(ones((s[0],1),M.dtype),v2 ** 2)),1)
    # inversion
    d=a11*a22 - a12 ** 2 + _esp
    ai11=a22 / d 
    ai12=- a12 / d
    ai22=a11 / d

    while ((np.linalg.norm(dp1) > tol) and (miter < opts.max_iter)):

        miter=miter + 1
        W= abs(M) * exp(1j*np.asarray((dot(p0,v1) + dot(p1,v2))))
        
        # projection
        pr1 = np.sum(np.conj(1j*W) * multiply(dot(ones((s[0],1),M.dtype),v1),(M-W)),axis=1)
        pr2 = np.sum(np.conj(1j*W) * multiply(dot(ones((s[0],1),M.dtype),v2),(M-W)),axis=1)
        
        dp0=np.real(ai11*pr1 + ai12*pr2).reshape([ai11.size,1])
        dp1=np.real(ai12*pr1 + ai22*pr2).reshape([ai11.size,1])
        dp1[np.isnan(dp1)]=0
        dp0[np.isnan(dp0)]=0 
        
        # update 
        p1=p1 + dp1
        p0=p0 + dp0

    
    dp1=sqrt(ai22)
    dp1[np.isnan(dp1)]=0
    dp1[np.isinf(dp1)]=0
    # relative residual
    res=M - multiply(abs(M),exp(1j*np.array(dot(p0,v1) + dot(p1,v2))))
    relres=np.sum(abs(res) ** 2,axis=1) / np.sum(abs(M) ** 2,axis=1)
    relres[np.isnan(relres)]=0
    
    # careful, mod two param should be the same data type
    # error when p1 is complex and pi is float
    p1[p1 > pi]=mod(p1[p1 > pi] + pi,2*pi) - pi
    p1[p1 < - pi]=mod(p1[p1 < - pi] + pi,2*pi) - pi
    
    p0=np.reshape(p0,s0[:(L_s0-1)])
    p1=np.reshape(p1,s0[:(L_s0-1)])
    dp1=np.reshape(dp1,s0[:(L_s0-1)])
    relres=np.reshape(relres,s0[:(L_s0-1)])
    
    return p1,dp1,relres,p0
