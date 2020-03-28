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

def QSMCalAll(mag, pha, voxel_size, B0_dir, TE, CF,isCSFCorr=True,Mask=[],qsm_filepath=''):
    
    iField =  mag * np.exp(-pha*1j)
    iField = iField.astype(np.complex64)
    
    matrix_size = np.array(iField.shape[:3])
#     voxel_size = np.array([0.9375,0.9375,2.0])
#     B0_dir = np.array([0,0,1])
#     CF = 123225803
#     TE=np.array([0.0036,0.0095,0.0154,0.0213,0.0273,0.0332,0.0391,0.0450])
    delta_TE = np.diff(TE).mean()

    print('run Fit_phase')
    iFreq_raw,N_std,a,b=Fit_ppm_complex(iField,10)
    print('run unwrapping')
    print(a)
    print(b)
    iFreq = unwrapLaplacian(iFreq_raw,matrix_size,voxel_size)
    
    ## do bet to remove brain skull using external exe 
    if Mask == []:
        qsm_maskpath = getBasepath(qsm_filepath) +'_mask.nii.gz'
        print('run Fit_ppm_complex')
        QSMBet(qsm_filepath, qsm_maskpath)
        
        fimg  = nib.load(qsm_maskpath)
        Mask = np.asarray(fimg.dataobj)[:,::-1,:]

    ######################################################################
    # Back ground field removal
    print('run background field removal')
    RDF, shim =PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir) # ,n_CG=2)
    
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

    outputVolume = [ iFreq_raw,  iFreq,  Mask_final     , RDF , QSM*1e6 , R2s]
    outputLabel  = ['iFreq_raw','iFreq','Mask_Brain+CSF','RDF','QSMx1e6','R2s']
    outputParam = dict(zip(outputLabel, outputVolume))
    return outputParam

def QSMBet(qsm_filepath,qsm_maskpath,f=0.3):
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