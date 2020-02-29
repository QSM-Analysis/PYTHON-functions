# Generated with SMOP  0.41
# import pandas as pd
import numpy as np
from dipole_kernel import dipole_kernel
from dipole_term import dipole_term
from cgsolve import cgsolve
# from sklearn.externals import joblib
# from dipole_kernel import*
# from cgsolve import*
# from dipole_term import*

# PDF.m

    # Projection onto Dipole Fields (PDF)
#   [RDF shim] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir,tol)
# 
#   output
#   RDF - the relative difference field, or local field
#   shim (optional) - the cropped background dipole distribution
# 
#   input
#   iFreq - the unwrapped field map
#   N_std - the noise standard deviation on the field map. (1 over SNR for single echo)
#   Mask - a binary 3D matrix denoting the Region Of Interest
#   matrix_size - the size of the 3D matrix
#   voxel_size - the size of the voxel in mm
#   B0_dir - the direction of the B0 field
#   tol(optional) - tolerance level
    
    #   When using the code, please cite 
#   T. Liu et al. NMR Biomed 2011;24(9):1129-36
#   de Rochefort et al. MRM 2010;63(1):194-206
    
    #   Created by Tian Liu in 2009
#   Modified by Tian Liu on 2011.02.01
#   Last modified by Tian Liu on 2013.07.24
    

def PDF(iFreq=None,N_std=None,Mask=None,matrix_size=None,voxel_size=None,B0_dir=None,tol=0.1,n_CG=30,space='imagespace',n_pad=40,*args,**kwargs):
    
# PDF.m:41
    
    # zero pad
    matrix_size0=np.copy(matrix_size)
    matrix_size0=matrix_size0.astype(int)
    matrix_size0=matrix_size0.flatten(1)
    #matrix_size0=matrix_size0.conj().T
# PDF.m:46
    d1=Mask.max(axis=1)
   
    #d1=Mask.reshape((int(matrix_size[0][0]),1,int(matrix_size[0][2])))
    d1=d1.max(axis=1)
    
# PDF.m:47
    d1first=(np.transpose(np.nonzero(d1)))[0]
# PDF.m:48
    d1last=(np.transpose(np.nonzero(d1)))[-1]
# PDF.m:49
    d2=((Mask.max(axis=2)).max(axis=0))
# PDF.m:51
    d2first=(np.transpose(np.nonzero(d2)))[0]
# PDF.m:52
    d2last=(np.transpose(np.nonzero(d2)))[-1]
# PDF.m:53
    d3=Mask.max(axis=0).max(axis=0)
# PDF.m:55
    d3first=(np.transpose(np.nonzero(d3)))[0]
# PDF.m:56
    d3last=((np.transpose(np.nonzero(d3)))[-1])
   
# PDF.m:57
    if n_pad > 0:
        matrix_size=np.concatenate([np.dot(np.floor((d1last - d1first + n_pad) / 2),2),np.dot(np.floor((d2last - d2first + n_pad) / 2),2),np.dot(np.floor((d3last - d3first + n_pad) / 2),2)])
# PDF.m:60
        #iFreq=iFreq[np.transpose(arange(d1first,d1last)),np.transpose(arange(d2first,d2last)),np.transpose(arange(d3first,d3last))] 
        iFreq=np.array(iFreq)[int(d1first[0]):int(d1last[0]),int(d2first[0]):int(d2last[0]),int(d3first[0]):int(d3last[0])]
# PDF.m:63
        N_std=N_std[int(d1first[0]):int(d1last[0]),int(d2first[0]):int(d2last[0]),int(d3first[0]):int(d3last[0])]
# PDF.m:64
        Mask=Mask[int(d1first[0]):int(d1last[0]),int(d2first[0]):int(d2last[0]),int(d3first[0]):int(d3last[0])]
# PDF.m:65
        padsize=([[0,int((matrix_size[0] - iFreq.shape[0]))],[0,int((matrix_size[1] - iFreq.shape[1]))],[0,int((matrix_size[2] - iFreq.shape[2]))]])
# PDF.m:66
        iFreq=np.pad(iFreq,padsize,'constant')
# PDF.m:67
        N_std=np.pad(N_std,padsize,'constant')
# PDF.m:68
        Mask=np.pad(Mask,padsize,'constant')
# PDF.m:69
    
    # generate the weighting
    W=1.0 / N_std
# PDF.m:74
    W[np.isinf(W)]=0
# PDF.m:75
    W=W*(Mask > 0)
    print('W',W.shape)
    print(W[100:102,100:102,50:52])
# PDF.m:76
    W_std=np.copy(W)
# PDF.m:77
    W_var=W ** 2
    print('W_var',W_var.shape)
    print(W_var[100:102,100:102,50:52])
# PDF.m:78
    ###### start the PDF method #####
    #if norm(ravel(B0_dir),1) < 1.01:
    

    if np.linalg.norm(B0_dir.flatten(1),1)<1.01:
        
        D=dipole_kernel([matrix_size,voxel_size,B0_dir])
# PDF.m:82
    else:
        D=dipole_kernel([matrix_size,voxel_size,B0_dir,space])
# PDF.m:84
    print('D',D.shape)
    print(D[100:102,100:102,50:52])
    # generating the RHS vector in Eq. 6 in the PDF paper
    p_temp=np.real(np.fft.ifftn(D*np.fft.fftn(W_var*(iFreq))))
    print('p_temp',p_temp.shape)
    print(p_temp[100:102,100:102,50:52])
# PDF.m:88
    b=p_temp[Mask == 0]
    print('b',b.shape)
    
# PDF.m:89
    # set erance level and maximum iteration allowed
    E_noise_level=np.real(np.fft.ifftn((D*np.fft.fftn((W_std*np.ones(np.shape(N_std)))))))
    print('E_noise_level',E_noise_level.shape)
    itermax=np.copy(n_CG)
    print('itermax=',itermax)
# PDF.m:93
    

    A=dipole_term
    joblib.dump(W,'W.pkl')
    joblib.dump(D,'D.pkl')
    joblib.dump(Mask,'Mask.pkl')

# PDF.m:94
    cg_tol=np.dot(tol,np.linalg.norm(E_noise_level[Mask == 0]/ np.linalg.norm(b.flatten(1))))
# PDF.m:95
    print ('cg_tol',cg_tol)
    x,res,num_iter=cgsolve(A,b,cg_tol,itermax,0)
# PDF.m:97
    print('CG stops at: res %f, iter %d\n',res,num_iter)
    print('x=',x)
    xx=np.zeros(D.shape)
    
# PDF.m:100
    xx[Mask == 0]=x
# PDF.m:101
    xx[Mask>0]=0
# PDF.m:102
    # background dipole field
    p_dipole=np.real(np.fft.ifftn((D*np.fft.fftn(xx))))
# PDF.m:105
    p_final=((iFreq - p_dipole)*Mask)
# PDF.m:108
    # remove zero pad
    if n_pad > 0:
        RDF=np.zeros(matrix_size0)
# PDF.m:112
        RDF[int(d1first):int(d1last),int(d2first):int(d2last),int(d3first):int(d3last)]=p_final[1:int(d1last- d1first + 1),1:int(d2last - d2first + 1),1:int(d3last - d3first + 1)]
# PDF.m:113
        shim=np.zeros(matrix_size0)
# PDF.m:116
        shim[int(d1first):int(d1last),int(d2first):int(d2last),int(d3first):int(d3last)]=xx[1:int(d1last - d1first + 1),1:int(d2last - d2first + 1),1:int(d3last - d3first + 1)]
# PDF.m:117
    else:
        RDF=np.copy(p_final)
# PDF.m:120
        shim=np.copy(xx)
# PDF.m:121
    return RDF
    