# Generated with SMOP  0.41
from smop.libsmop import *
#import tensorflow as tf
import sys
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
    
import numpy as np
from numpy import zeros, min, pi, abs,ones,dot, multiply,sum

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
        sz=shape(M)
        numpy.zeros
        p1=zeros(M.shape[:3],M.dtype)
        dp1=zeros(M.shape[:3],M.dtype)
        relres=zeros(M.shape[:3],M.dtype)
        p0=zeros(M.shape[:3],M.dtype)
        n=ceil(bytes / 1000000000.0)
        ns=floor(sz(3) / n)
        for slice in arange(1,sz(3),ns).reshape(-1):
            rng=arange(slice,min(slice + ns - 1,sz(3)))
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
            M=sum(multiply(M,conj(repmat(M[:,:,:,1,:],concat([1,1,1,shape(M,4),1])))),5)
            M=multiply(sqrt(abs(M)),exp(dot(1j,angle(M))))
    
    # determine angle
    M=np.conjugate(M)
    s0=M.shape
    L_s0=len(s0)
    if L_s0<=len(M.shape):
        nechos=M.shape[L_s0-1]
    else:
        nechos=1
    # M=np.reshape(M,[np.prod(M.shape[:-1]),M.shape[-1]])
    M=np.reshape(M,[np.prod(s0[:(L_s0-1)]),s0[L_s0-1]])
    s=M.shape
    Y=np.angle(M[:,:np.min([3,nechos])])
    c=Y[:,1] - Y[:,0]
    temp = np.array([abs(c-2*pi),abs(c),abs(c+2*pi)])
    m = temp.min(1)
    ind = temp.argmin(1)
    k=ind
    for i in range(len(ind)):
        if ind[i]==0:
            k[i]=1
        else:
            k[i]=0
            
    c[k] = c[k] - 2*pi
    
    for i in range( len(ind)):
        if ind[i]==2:
            k[i]=1
        else:
            k[i]=0
    c[k] = c[k] + 2*pi
    
    for n in range(min([2,nechos-1])):
        cd=((Y[:,n + 1] - Y[:,n])) - c
        
        
        Y[cd < - pi,(n + 1):len(Y[cd<-pi])]=Y[cd < - pi,(n + 1):len(Y[cd<-pi])] + 2*pi
        Y[cd> pi,(n + 1):len(Y[cd>pi])]=Y[cd> pi,(n + 1):len(Y[cd>pi])] + 2*pi
       #Y[cd < - pi,arange((n + 1),len(Y[cd<-pi]))]=Y[cd < - pi,arange((n + 1),len(Y[cd<-pi]))] + 2*pi
       # for i in range((n + 1),len(Y[cd<-pi])):
           # Y[cd < - pi,i]+=dot(2,pi)
# Fit_ppm_complex.m:92
        #for i in range((n + 1),len(Y[cd>pi])): 
            #Y[cd > pi,i]=Y[cd > pi,i] - dot(2,pi)
# Fit_ppm_complex.m:93
  #A=tf.cast(concat([[1,0],[1,1],[1,2]]),M.dtype)
    A=np.array([[1,0],[1,1],[1,2]])
# Fit_ppm_complex.m:95
    #ip = A(1:min(3,nechos),:)\Y(:,1:min(3,nechos))';
    ip=dot(np.mat(A[:min([3,nechos]),:]).I,Y[:,:min([3,nechos])].conj().T)

# Fit_ppm_complex.m:96
    p0=np.array(ip[0,:]).conj().T
# Fit_ppm_complex.m:97
    p1=np.array(ip[1,:]).conj().T
# Fit_ppm_complex.m:98
    dp1=copy(p1)
# Fit_ppm_complex.m:100
    tol=dot(np.linalg.norm(p1[:]),opts.reltol)
# Fit_ppm_complex.m:101
    iter=0
# Fit_ppm_complex.m:102
    # max_iter = 30;
    
    # weigthed least square
# calculation of WA'*WA
    numpy.ones
    
    v1=ones((1,nechos),M.dtype)
    
# Fit_ppm_complex.m:107
    #v2=cast((arange(0,(nechos - 1))),M.dtype)
    v2=np.array([[i for i in range(nechos)]])

# Fit_ppm_complex.m:108
    a11=sum(multiply(abs(M) ** 2.0,(dot(ones((s[0],1),M.dtype),(v1 ** 2)))),2)
# Fit_ppm_complex.m:114
    a12=sum(multiply(abs(M) ** 2.0,(dot(ones((s[0],1),M.dtype),(multiply(v1,v2))))),2)
# Fit_ppm_complex.m:115
    a22=sum(multiply(abs(M) ** 2.0,(dot(ones((s[0],1),M.dtype),(v2 ** 2)))),2)
# Fit_ppm_complex.m:111
    # inversion
    d=multiply(a11,a22) - a12 ** 2
# Fit_ppm_complex.m:113
    ai11=a22 / d 
# Fit_ppm_complex.m:114
    ai12=- a12 / d
# Fit_ppm_complex.m:115
    ai22=a11 / d
# Fit_ppm_complex.m:116
    while ((np.linalg.norm(dp1) > tol) and (iter < opts.max_iter)):

        iter=iter + 1
# Fit_ppm_complex.m:119
        W=abs(M)*exp(multiply(1j,(dot(p0,v1) + dot(p1,v2))))
        
# Fit_ppm_complex.m:120
        pr1=sum(multiply(multiply(np.conj(dot(1j,W)),(dot(ones((s[0],1),M.dtype),v1))),(M - W)),2)
        #pr1=sum(np.conj(dot(1j,W))*(multiply(ones((s[0],1),M.dtype),v1))*(M - W),2)
        
# Fit_ppm_complex.m:123
        pr2=sum(multiply(multiply(np.conj(dot(1j,W)),(dot(ones((s[0],1),M.dtype),v2))),(M - W)),2)
# Fit_ppm_complex.m:124
        dp0=np.array([np.real(multiply(ai11,pr1) + multiply(ai12,pr2))]).T
# Fit_ppm_complex.m:126
        dp1=np.array([np.real(multiply(ai12,pr1) + multiply(ai22,pr2))]).T
# Fit_ppm_complex.m:127
        dp1[np.isnan(dp1)]=0
# Fit_ppm_complex.m:128
        dp0[np.isnan(dp0)]=0 
# Fit_ppm_complex.m:129
        p1=p1 + dp1
# Fit_ppm_complex.m:132
        p0=p0 + dp0
# Fit_ppm_complex.m:133

    
    # error propagation
    dp1=sqrt(ai22)
# Fit_ppm_complex.m:139
    dp1[np.isnan(dp1)]=0
# Fit_ppm_complex.m:140
    dp1[np.isinf(dp1)]=0
# Fit_ppm_complex.m:141
    # relative residual
    res=M - multiply(abs(M),exp(dot(1j,(dot(p0,v1) + dot(p1,v2)))))
# Fit_ppm_complex.m:144
    relres=sum(abs(res) ** 2,2) / sum(abs(M) ** 2,2)
# Fit_ppm_complex.m:145
    relres[np.isnan(relres)]=0
# Fit_ppm_complex.m:146
    p1[p1 > pi]=mod(p1[p1 > pi] + pi,dot(2,pi)) - pi
# Fit_ppm_complex.m:148
    p1[p1 < - pi]=mod(p1[p1 < - pi] + pi,dot(2,pi)) - pi
# Fit_ppm_complex.m:149
    p0=np.reshape(p0,s0[:(L_s0-1)])
# Fit_ppm_complex.m:151
    p1=np.reshape(p1,s0[:(L_s0-1)])
# Fit_ppm_complex.m:152
    dp1=np.reshape(dp1,s0[:(L_s0-1)])
# Fit_ppm_complex.m:153
    relres=np.reshape(relres,s0[:(L_s0-1)])
    return p1,dp1,relres,p0
# Fit_ppm_complex.m:154