# -*- coding: utf-8 -*-
#from smop.libsmop import *
# extract_CSF.m
import numpy as np
import scipy.io as sio
import scipy.ndimage as ndi
from skimage import measure,color
from numpy import ravel, sum, sqrt, abs, multiply
from scipy.ndimage.morphology import binary_erosion
from skimage.morphology.misc import remove_small_objects
from PyQSM.SMV import SMV
#from SMV import SMV

def isempty(data):
    if (type(data)!=np.ndarray):
        if (data == None) or (data == []):
            return True
    return False

def extract_CSF(R2s=None,Mask=None,voxel_size=None,flag_erode=1,thresh_R2s=5,*args,**kwargs):
#     varargin = extract_CSF.varargin
#     nargin = extract_CSF.nargin

    if isempty(R2s):
        Mask_ROI_CSF=[]
        return Mask_ROI_CSF
    
    n_region_cen=3
    matrix_size=np.array(Mask.shape)
    
    X,Y,Z = np.mgrid[0:matrix_size[0],0:matrix_size[1],0:matrix_size[2]]
    X, Y, Z = X * voxel_size[0], Y * voxel_size[1], Z * voxel_size[2]
    X_cen=sum(X[Mask > 0]) / sum(ravel(Mask))
    Y_cen=sum(Y[Mask > 0]) / sum(ravel(Mask))
    Z_cen=sum(Z[Mask > 0]) / sum(ravel(Mask))

    # Center region (sphere)
    radius_cen=30
    Mask_cen=sqrt(abs(X - X_cen) ** 2 + abs(Y - Y_cen) ** 2 + abs(Z - Z_cen) ** 2) <= radius_cen
    if flag_erode:
        # temporally use erosion to replace SMV
        #stru = np.ones((np.array(Mask.shape)*0.05).astype(int))
        #Mask = binary_erosion(Mask,structure=stru.astype(Mask.dtype))
        Mask=SMV(Mask,[matrix_size,voxel_size,10]) > 0.999

    # find csf condidate region in mask center circle    
    Mask_raw_1=(R2s < thresh_R2s) * Mask_cen
    # additional open to remove small region
    Mask_raw_1 = remove_small_objects(Mask_raw_1,5)#,np.ones([1,1,1]).astype(Mask_raw_1.dtype))
    def getCCind(Mask):
        CC,num=measure.label(Mask,return_num=True,neighbors=4)  
        CCsize= []
        for i in range(num):
            CCsize.append(np.sum(CC==i))
        CCind = np.argsort(CCsize)[::-1] # descending sort
        return CC, CCind
    CC, CCind = getCCind(Mask_raw_1)
    ROIs_region_cen = np.zeros(matrix_size)
    if CCind.size < n_region_cen:
        return [] 
    for i in range(1,n_region_cen+1):# not include background i ==0
        ROIs_region_cen[CC==CCind[i]] = i 


    # find csf candidate region in mask   
    Mask_raw_2=(R2s < thresh_R2s) * Mask
    # additional open to remove small region
    Mask_raw_2 = remove_small_objects(Mask_raw_2,5)#opening(Mask_raw_2)#,np.ones([3,3,3]).astype(Mask_raw_2.dtype))
    
    CC, CCind = getCCind(Mask_raw_2)
    ROIs_region = np.zeros(matrix_size)
    for i in range(1,CCind.size):# not include background i ==0
        ROIs_region[CC==CCind[i]] = i 
        
    # Choose regions which appear at center
    Mask_ROI_CSF=np.zeros(matrix_size)
    for i in np.unique(ROIs_region[(ROIs_region_cen >0) * (ROIs_region > 0)]):
        Mask_ROI_CSF[ROIs_region == i]=1
    
    Mask_ROI_CSF[Mask == 0]=0

    return Mask_ROI_CSF
    
if __name__ == '__main__':
    pass
# import matplotlib.pyplot as plt
# fig=plt.figure()
# fig1=fig.add_subplot(231)
# fig1.imshow(Mask[:,:,30],'gray',vmin=0,vmax=5)

















