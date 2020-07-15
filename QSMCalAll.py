# -*- coding: utf-8 -*-
# Created on Mar 10, 2020
# Author: Maxwell4444
# All right reserved

import sys
import numpy as np
from PyQSM.Fit_ppm_complex import Fit_ppm_complex
from PyQSM.unwrapLaplacian import unwrapLaplacian
import os
import nibabel as nib
from PyQSM.arlo import arlo
from PyQSM.extract_CSF import extract_CSF, isempty
from PyQSM.PDF import PDF
from PyQSM.MEDI_L1 import MEDI_L1
from Apps.MyTools.Common import getBasepath
from MPApp.DataReaderPackage.DataReader import DataReader
from MPApp.Seg.SegAlg import BinaryErode

def QSMCalAll(mag, pha, voxel_size, B0_dir, TE, CF,isCSFCorr=True,Mask=[],qsm_filepath=''):
    
    mag = mag.swapaxes(0,1)[:,:,::-1]
    pha = pha.swapaxes(0,1)[:,:,::-1]
    Mask = Mask.swapaxes(0,1)[:,:,::-1]
    voxel_size = voxel_size[[1,0,2]]
    #B0_dir = B0_dir[[1,0,2]]
    
    iField =  mag * np.exp(-pha*1j)
    iField = iField.astype(np.complex64)
    
    matrix_size = np.array(iField.shape[:3])
#     voxel_size = np.array([0.9375,0.9375,2.0])
#     B0_dir = np.array([0,0,1])
#     CF = 123225803
#     TE=np.array([0.0036,0.0095,0.0154,0.0213,0.0273,0.0332,0.0391,0.0450])
    delta_TE = np.diff(TE).mean()

    ################################################################################
    # import MEDI matlab results for comparison
#     import scipy.io as scio
#     data = scio.loadmat('D:\work\MATLAB\!QSM\code\MEDI_toolbox_201904\qsm_data.mat')
#     iField_mat = data['iField'].astype(np.complex64)
#     iFreq_raw_mat = data['iFreq_raw'].astype(float)#.swapaxes(0,1)[:,:,::-1]
#     N_std_mat = data['N_std'].astype(float)#.swapaxes(0,1)[:,:,::-1]
#     iFreq_mat = data['iFreq'].astype(float)#.swapaxes(0,1)[:,:,::-1]
#     Mask_mat = data['Mask'].astype(float)#.swapaxes(0,1)[:,:,::-1]
#     RDF_mat = data['RDF'].astype(float)#.swapaxes(0,1)[:,:,::-1]
#     Mask_CSF_mat = data['Mask_CSF'].astype(float)
# #     B0_dir_mat = data['B0_dir'].astype(float)
    ################################################################################
    
    print('run Fit_phase')
    iFreq_raw,N_std,a,b=Fit_ppm_complex(iField)#,5)
#     iFreq_raw = iFreq_raw_mat
#     N_std = N_std_mat
    
    print('run unwrapping')
    iFreq = unwrapLaplacian(iFreq_raw,matrix_size,voxel_size)

    ################################################################################
    # display iFreq_raw, iFreq and compared with MEDI matlab code
#     import matplotlib.pyplot as plt
#     fig=plt.figure()
#     fig1=fig.add_subplot(231)
# #     fig1.imshow(iFreq_raw[:,:,30],'gray',vmin=-1,vmax=1)
#     fig2=fig.add_subplot(232)
# #     fig2.imshow(iFreq[:,:,30],'gray',vmin=-1,vmax=1)
#     fig1.imshow(np.abs(iFreq_raw_mat-iFreq_raw)[:,:,40],'gray',vmin=-1,vmax=1)
#     fig2.imshow(np.abs(iFreq_mat-iFreq)[:,:,40],'gray',vmin=-1,vmax=1)
    ################################################################################
    
    ## do bet to remove brain skull using external exe 
    if Mask == []:
        qsm_maskpath = getBasepath(qsm_filepath) +'_mask.nii.gz'
        print('run Fit_ppm_complex')
        QSMBet(qsm_filepath, qsm_maskpath,f=0.5)
        
#         fimg  = nib.load(qsm_maskpath)
#         Mask = np.asarray(fimg.dataobj)[:,::-1,:]
        datareader = DataReader()
        datareader.LoadNiftyFile(qsm_maskpath)
        Mask = datareader._data
    #Mask = BinaryErode(Mask,structure=np.ones([2]*3))
    
    # Back ground field removal)
    print('run background field removal')
    RDF, shim =PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir) # ,n_CG=2)
    #RDF_matcal, shim =PDF(iFreq.swapaxes(0,1)[:,:,::-1],N_std.swapaxes(0,1)[:,:,::-1],Mask.swapaxes(0,1)[:,:,::-1],matrix_size[[1,0,2]],voxel_size[[1,0,2]],data['B0_dir'].astype(float))
    
    ######################################################################
    # display RDF and Mask
#     fig3=fig.add_subplot(233)
#     fig3.imshow(Mask[:,:,30],'gray',vmin=0,vmax=1)
#     fig4=fig.add_subplot(234)
#     fig4.imshow(RDF[:,:,30],'gray',vmin=-1,vmax=1)
#      
#     fig4.imshow(np.abs(RDF_mat-RDF)[:,:,40],'gray',vmin=-1,vmax=1)
#     fig3.imshow(np.abs(Mask_mat-Mask)[:,:,40],'gray',vmin=-1,vmax=1)
    
    if isCSFCorr:
        # R2s cal
        print('run R2s')
        R2s=arlo(TE/1000,mag)*Mask
        #R2s=arlo(TE,abs(iField))
        # CSF Mask
        print('run Mask_CSF')
        Mask_CSF=extract_CSF(R2s,Mask,voxel_size)
        if not isempty(Mask_CSF):
            Mask_CSF = Mask_CSF>0
    else:
        Mask_CSF = []
        R2s = []
    # cal MEDI_L1
    # iMag=np.sqrt(np.sum(abs(iField) ** 2,3))
    iMag=np.sqrt(np.sum(mag ** 2,3))
    print('run Dipole inversion')
    QSM,cost_reg_history,cost_data_history = MEDI_L1(RDF,N_std,iMag,Mask,Mask_CSF, voxel_size,B0_dir,CF, delta_TE, _lambda=1000, edge_percentage=0.9,smv_radius=5)#,max_iter=1,cg_max_iter=1)
    #QSM = iMag
    
    if isempty(Mask_CSF): Mask_final = Mask
    else: Mask_final = Mask+Mask_CSF

    outputVolume = [ iFreq_raw,  iFreq,  Mask_final     , RDF , -QSM*1e6 , R2s]
    for i in range(len(outputVolume)):
        outputVolume[i] = outputVolume[i].swapaxes(0,1)[:,:,::-1]
    outputLabel  = ['iFreq_raw','iFreq','Mask+CSF','RDF','QSM','R2s']
    outputParam = dict(zip(outputLabel, outputVolume))
    return outputParam

def QSMBet(qsm_filepath,qsm_maskpath,f=0.5):
    options=' -m ' + qsm_maskpath + ' -f ' + str(f) + ' -g 0'
    cmd = ''
    for basepath in ['..\\PyQSM\\','']:
        curpath = basepath + 'bet2.exe'
        if os.path.isfile(curpath):
            cmd = curpath +  ' ' + qsm_filepath + ' ' + getBasepath(qsm_filepath)+'_masked.nii.gz' + options
    if cmd == '': 
        print('bet2.exe not found')
        return
    os.system(cmd)
    
def getBasepath(f):
    basepath = f
    if f[-4:] == '.npz': basepath = f[:-4]
    elif f[-4:] == '.nii': basepath = f[:-4]
    elif f[-7:] == '.nii.gz': basepath = f[:-7]
    else:basepath = f.split('.')[0]
    return basepath

if __name__ == '__main__':
    qsm_filepath = r''
    # load data
    fimg  = nib.load(qsm_filepath)
    data = np.asarray(fimg.dataobj)