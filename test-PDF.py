# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:19:16 2020

@author: TR
"""
from PDF import*

iFreq=scio.loadmat('matlabdata_iFreq_FOR RDF.mat')['iFreq']
Mask=scio.loadmat('matlabdata_Mask_FOR RDF.mat')['Mask']
matrix_size=scio.loadmat('matlabdata_matrix_size_FOR RDF.mat')['matrix_size']
voxel_size=scio.loadmat('matlabdata_voxel_size_FOR RDF.mat')['voxel_size']
B0_dir=scio.loadmat('matlabdata_ B0_dir_FOR RDF.mat')['B0_dir']
N_std=scio.loadmat('matlabdata_N_std_FOR RDF.mat')['N_std']

RDF=PDF(iFreq,N_std,Mask,matrix_size,voxel_size,B0_dir)