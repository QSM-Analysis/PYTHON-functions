# Generated with SMOP  0.41
#from smop.libsmop import *
# extract_CSF.m
import numpy as np
import scipy.io as sio
from numpy import*
from SMV import*
import scipy.ndimage as ndi
from skimage import measure,color

def isempty(data):
    if (type(data)!=np.ndarray):
        if (data == None) or (data == []):
            return True
    return False

#@function
def extract_CSF(R2s=None,Mask=None,voxel_size=None,flag_erode=1,thresh_R2s=5,*args,**kwargs):
#     varargin = extract_CSF.varargin
#     nargin = extract_CSF.nargin

    if isempty(R2s):
        Mask_ROI_CSF=[]
        return Mask_ROI_CSF
    
    n_region_cen=3
    matrix_size=np.array(np.shape(Mask))
    
    X,Y,Z = np.mgrid[0:matrix_size[0],0:matrix_size[1],0:matrix_size[2]]
    X, Y, Z = X * voxel_size[0], Y * voxel_size[1], Z * voxel_size[2]
    X_cen=sum(X[Mask > 0]) / sum(ravel(Mask))
    Y_cen=sum(Y[Mask > 0]) / sum(ravel(Mask))
    Z_cen=sum(Z[Mask > 0]) / sum(ravel(Mask))

    # Center region (sphere)
    radius_cen=30
# extract_CSF.m:27
    Mask_cen=sqrt(abs(X - X_cen) ** 2 + abs(Y - Y_cen) ** 2 + abs(Z - Z_cen) ** 2) <= radius_cen
# extract_CSF.m:28
    if flag_erode:
        Mask=SMV(Mask,[[matrix_size],voxel_size,10]) > 0.999
# extract_CSF.m:32
    
    Mask_raw_1=multiply((R2s < thresh_R2s),Mask_cen)
    print(Mask_raw_1.shape)
# extract_CSF.m:37

#    CC=bwconncomp(Mask_raw_1,6)      #获得6连通区域      smop
    CC,num=measure.label(Mask_raw_1,return_num=True,connectivity=2)  #如果是二维图像，则是8连通区域  STILL HAVE PROBLEMS
# extract_CSF.m:38

#    numPixels=cellfun(numel,CC.PixelIdxList)    #计算联通区域面积      
    numPixels={}
    for i in range(1,num+1):
        numPixels[i+1]=np.size(np.where(CC==i)[1])  
# extract_CSF.m:39

    numPixels_sorted=sorted(numPixels,reverse = True)        #对各行元素降序排列
    idxs=np.argsort(numPixels.values,axis=0)                      #使numPixels_sorted[idxs]=numPixels
# extract_CSF.m:40
   
    ROIs_region_cen=np.zeros(matrix_size.shape)
    CC_PixelIdxList=np.arange(CC.size)
# extract_CSF.m:41
    for i in np.arange(1,n_region_cen).reshape(-1):
        numPixels_sorted[i]
        idx=idxs[i]
# extract_CSF.m:44     
        ROIs_region_cen[CC_PixelIdxList[idx]]=i
# extract_CSF.m:45
    
    Mask_raw_2=multiply((R2s < thresh_R2s),Mask)
# extract_CSF.m:48
    CC,num=measure.label(Mask_raw_2,return_num=True,connectivity=2)
# extract_CSF.m:49
    numPixels={}
    for i in range(1,num+1):
        numPixels[i+1]=np.size(np.where(CC==i)[1])
# extract_CSF.m:50
    numPixels_sorted=sorted(numPixels,reverse = True)
    idxs=np.argsort(numPixels.values,axis=0) 
# extract_CSF.m:51
    ROIs_region=np.zeros(matrix_size.shape)
# extract_CSF.m:52
    for i in arange(1,np.length(idxs)).reshape(-1):
        numPixels_sorted(i)
        idx=idxs(i)
# extract_CSF.m:55
        ROIs_region[CC.PixelIdxList[idx]]=i
# extract_CSF.m:56
    
    # Choose regions which appear at center
    Mask_ROI_CSF=zeros(matrix_size)
# extract_CSF.m:60
    for i in transpose(unique(ROIs_region(ROIs_region_cen > logical_and(0,ROIs_region) > 0))).reshape(-1):
        Mask_ROI_CSF[ROIs_region == i]=1
# extract_CSF.m:62
    
    Mask_ROI_CSF[Mask == 0]=0
# extract_CSF.m:65
    return Mask_ROI_CSF
    
if __name__ == '__main__':
    pass


Mask=sio.loadmat('x00000001.mat')['Mask']
voxel_size=sio.loadmat('voxel_size.mat')['voxel_size']
R2s=sio.loadmat('R2s.mat')['R2s']
Mask_ROI_CSF = extract_CSF(R2s, Mask, voxel_size)
















