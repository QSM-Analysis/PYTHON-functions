# Generated with SMOP  0.41
import sys
from smop.libsmop import *
from Fit_ppm_complex import *
from unwrapLaplacian import *
# from BET import*
from arlo import*
# from extract_CSF import*
# from PDF import*
# from sphere_kernel import*
# from dipole_term import*
# from cgsolve import*
# from store_QSM_results import*
# from cgsolve import*
# from gradient_mask import*
# from cgrad import*
# from cdiv import*
# from bdiv import*
# from fgrad import*
# from dataterm_mask import*
# from SMV import*
# from sphere_kernel import*
# from SMV_kernel import*
import scipy.io as scio
import matplotlib.pyplot as plt
import numpy as np
# MEDI_Batch_yx.m

   
   # loadmat(medi_siemens_data)
    #print(iField)
data=scio.loadmat('medi_siemens_data.mat')
iField=data['iField']  

# iMag=sqrt(sum(abs(iField) ** 2,4))

# iFreq_raw,N_std,a,b=Fit_ppm_complex(iField)
# fig=plt.figure()
# fig1=fig.add_subplot(111)
# fig1.imshow(iFreq_raw[:,:,30])
# plt.show()
# data_iFreq_raw=scio.savemat('data_iFreq_raw.mat',{'iFreq_raw':iFreq_raw})
# 
#
data_iFreq_raw=scio.loadmat('matlabdata_iFreq_raw.mat')
data_iFreq=scio.loadmat('matlabdata_iFreq.mat')

iFreq_raw=data_iFreq_raw['iFreq_raw']
matrix_size=data['matrix_size'][0]
voxel_size=np.array((data['voxel_size']).T[0])

iFreq=unwrapLaplacian(iFreq_raw,matrix_size,voxel_size)
print(iFreq[128:132,128:132,30])
print(data_iFreq['iFreq'][128:132,128:132,30])

# fig2=fig.add_subplot(212)
# fig2.imshow(np.real(iFreq[:,:,30]),'gray')
# plt.show()
 
# Mask=BET(iMag,matrix_size,voxel_size)
# qsm_filepath='D:\eclipse-workspace\OncoImageAnalysis\src\QSM\data\004_QSM_MONO_8TE_IPAT2_68phpf_NoFC_mask.nii'

# qsm_filepath_mask='D:\eclipse-workspace\OncoImageAnalysis\src\QSM\data\004_QSM_MONO_8TE_IPAT2_68phpf_NoFC_matmask.nii'

# TE=data['TE']
# R2s=arlo(TE,abs(iField))

# Mask_CSF=extract_CSF(R2s,Mask,voxel_size)

# RDF=PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir)

# figure
# imshow(RDF(arange(),arange(),30),concat([- 1,1]))
# QSM=MEDI_L1('lambda',1000,'percentage',0.9,'smv',5)
# figure
# imshow(QSM(arange(),arange(),30),[])