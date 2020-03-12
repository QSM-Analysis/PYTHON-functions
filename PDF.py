# Generated with SMOP  0.41
# import pandas as pd
import numpy as np
from PyQSM.dipole_kernel import dipole_kernel
from PyQSM.dipole_term import dipole_term
from PyQSM.cgsolve import cgsolve
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
    
    # zero pad
    matrix_size0=matrix_size.astype(int)
    d1=Mask.max(axis=2).max(axis=1)
    d1first, d1last = np.nonzero(d1)[0][0], np.nonzero(d1)[0][-1]

    d2=Mask.max(axis=2).max(axis=0)
    d2first, d2last = np.nonzero(d2)[0][0], np.nonzero(d2)[0][-1]

    d3=Mask.max(axis=1).max(axis=0)
    d3first, d3last = np.nonzero(d3)[0][0], np.nonzero(d3)[0][-1]
    
    d1last1, d2last1,d3last1 = d1last+1, d2last+1, d3last+1
    if n_pad > 0:
        matrix_size=np.array([np.floor((d1last - d1first + n_pad) / 2)*2,
                                    np.floor((d2last - d2first + n_pad) / 2)*2,
                                    np.floor((d3last - d3first + n_pad) / 2)*2]).astype(int)
        
        iFreq = iFreq[d1first:d1last1, d2first:d2last1, d3first:d3last1]
        N_std = N_std[d1first:d1last1, d2first:d2last1, d3first:d3last1]
        Mask  =  Mask[d1first:d1last1, d2first:d2last1, d3first:d3last1]
        
        padsize=([[0, matrix_size[0] - iFreq.shape[0]],
                  [0, matrix_size[1] - iFreq.shape[1]],
                  [0, matrix_size[2] - iFreq.shape[2]]])
        iFreq = np.pad(iFreq,padsize,'constant')
        N_std = np.pad(N_std,padsize,'constant')
        Mask  = np.pad(Mask, padsize,'constant')
    
    # generate the weighting
    W=1.0 / N_std
    W[np.isinf(W)]=0
    W=W*(Mask > 0)
    W_std=W
    W_var=W ** 2

    ###### start the PDF method #####
    #if norm(ravel(B0_dir),1) < 1.01:
    if np.linalg.norm(B0_dir.flatten(1),1)<1.01:
        D=dipole_kernel(matrix_size,voxel_size,B0_dir)
    else:
        D=dipole_kernel(matrix_size,voxel_size,B0_dir,space)

    # generating the RHS vector in Eq. 6 in the PDF paper
    p_temp=np.real(np.fft.ifftn(D*np.fft.fftn(W_var*(iFreq))))

    b = p_temp[Mask == 0] #p_temp.flatten(1)[Mask.flatten(1)==0]#
    
    # set erance level and maximum iteration allowed
    E_noise_level=np.real(np.fft.ifftn((D*np.fft.fftn((W_std*np.ones(np.shape(N_std)))))))
    itermax=np.copy(n_CG)
    print('itermax=',itermax)

    #A=(dipole_term, W_var,D,Mask)
    A = lambda x=None: dipole_term(W_var, D,Mask, x)
    cg_tol=tol*np.linalg.norm(E_noise_level[Mask == 0]/ np.linalg.norm(b))
    print ('cg_tol',cg_tol)
    
    x,res,num_iter=cgsolve(A,b,cg_tol,itermax,0)
    print('CG stops at: res '+str(res)+' , iter ' + str(num_iter),'\n')
    #print('x=',x)
    
    xx=np.zeros(D.shape)
    xx[Mask == 0]=x
    xx[Mask>0]=0

    # background dipole field
    p_dipole=np.real(np.fft.ifftn((D*np.fft.fftn(xx))))
    p_final=((iFreq - p_dipole)*Mask)

    # remove zero pad
    if n_pad > 0:
        RDF=np.zeros(matrix_size0)
        RDF[int(d1first):int(d1last1),int(d2first):int(d2last1),int(d3first):int(d3last1)]=p_final[:int(d1last1- d1first),:int(d2last1 - d2first),:int(d3last1 - d3first)]
        shim=np.zeros(matrix_size0)
        shim[int(d1first):int(d1last1),int(d2first):int(d2last1),int(d3first):int(d3last1)]=xx[:int(d1last1 - d1first),:int(d2last1 - d2first),:int(d3last1 - d3first)]
    else:
        RDF=p_final
        shim=xx
    return RDF, shim  
    
# import matplotlib.pyplot as plt
# fig=plt.figure()
# fig1=fig.add_subplot(221)
# fig1.imshow(p_temp[:,:,50],'gray',vmin=0,vmax=1000000)
