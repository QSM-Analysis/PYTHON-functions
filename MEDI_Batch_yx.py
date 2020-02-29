import sys
import scipy.io as scio
import matplotlib.pyplot as plt
import numpy as np
from Fit_ppm_complex import Fit_ppm_complex
from unwrapLaplacian import unwrapLaplacian
import os
import nibabel as nib
from arlo import arlo
from PyQSM.extract_CSF import extract_CSF
from PyQSM.PDF import PDF

fig=plt.figure()

data=scio.loadmat('medi_siemens_data.mat')
iField=data['iField'].astype(np.complex64)

matrix_size = iField.shape[:3]
voxel_size = [0.9375,0.9375,2.0]
    
#iMag=np.sqrt(np.sum(abs(iField) ** 2,3))

######################################################################
run_Fit_ppm_complex=False
if run_Fit_ppm_complex:
    iFreq_raw,N_std,a,b=Fit_ppm_complex(iField)
    np.savez(r'QSM_temp_data', fit_ppm_complex=(iFreq_raw, N_std,a,b))
else:
    arch = np.load(r'QSM_temp_data.npz')
    iFreq_raw,N_std,a,b = arch['fit_ppm_complex']
fig1=fig.add_subplot(231)
fig1.imshow(abs(iFreq_raw[:,:,30]),'gray',vmin=0,vmax=1)

######################################################################
run_unwrapLaplacian = False
if run_unwrapLaplacian:
    
    
    iFreq = unwrapLaplacian(iFreq_raw,matrix_size,voxel_size)
    np.savez(r'QSM_temp_data2', unwraplaplacian=(iFreq))
    print(iFreq[128:132,128:132,30])
    # matlab_iFreq=scio.loadmat('matlabdata_iFreq.mat')
    # print(matlab_iFreq['iFreq'][128:132,128:132,30])
else:
    arch = np.load(r'QSM_temp_data2.npz')
    iFreq = arch['unwraplaplacian']
fig2=fig.add_subplot(232)
fig2.imshow(np.abs(iFreq[:,:,30]),'gray',vmin=0,vmax=1)

######################################################################
## bet 
## 1) do bet using external exe
# # qsm raw data 
# qsm_filepath = r'data\004_QSM_MONO_8TE_IPAT2_68phpf_NoFC.nii.gz'
# # bet, remove brain skull
# qsm_maskpath = qsm_filepath[:-7] +'_mask.nii.gz'
# options=' -m ' + qsm_maskpath + ' -f .5 -g 0'
# cmd = r'D:\work\MATLAB\!QSM\code\MEDI_toolbox_201904\bet2.exe ' + qsm_filepath + ' ' + qsm_filepath[:-7]+'_masked.nii.gz' + options
# os.system(cmd)

## 2) load bet data
maskpath = r'004_QSM_MONO_8TE_IPAT2_68phpf_NoFC_mask.nii.gz'
fimg  = nib.load(maskpath)
Mask = np.asarray(fimg.dataobj)[:,::-1,:]
# fig3=fig.add_subplot(233)
# fig3.imshow(Mask[:,:,30],'gray',vmin=0,vmax=1)

######################################################################
# Back ground field removal
B0_dir = np.array([0,0,1])
RDF=PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir)
fig3=fig.add_subplot(233)
fig3.imshow(RDF[:,:,30],'gray',vmin=-1,vmax=1)

######################################################################
# R2s cal
TE=np.array([0.0036,0.0095,0.0154,0.0213,0.0273,0.0332,0.0391,0.0450])
R2s=arlo(TE,abs(iField))
fig4=fig.add_subplot(234)
fig4.imshow(R2s[:,:,30]*Mask[:,:,30],'gray',vmin=0,vmax=40)

Mask_CSF=extract_CSF(R2s,Mask,voxel_size)
fig5=fig.add_subplot(235)
fig5.imshow(Mask_CSF[:,:,30],'gray',vmin=0,vmax=1)

# RDF=PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir)

# figure
# imshow(RDF(arange(),arange(),30),concat([- 1,1]))
# QSM=MEDI_L1('lambda',1000,'percentage',0.9,'smv',5)
# figure
# imshow(QSM(arange(),arange(),30),[])
plt.show()
pass