# Generated with SMOP  0.41
from smop.libsmop import *
# extract_CSF.m
import numpy as np
from numpy import*
    
@function
def extract_CSF(R2s=None,Mask=None,voxel_size=None,flag_erode=None,thresh_R2s=None,*args,**kwargs):
    varargin = extract_CSF.varargin
    nargin = extract_CSF.nargin

    if isempty(R2s):
        Mask_ROI_CSF=[]
# extract_CSF.m:4
        return Mask_ROI_CSF
  
    if nargin < 5:
        thresh_R2s=5
# extract_CSF.m:9
    
    if nargin < 4:
        flag_erode=1
# extract_CSF.m:12
    
    n_region_cen=3
# extract_CSF.m:15
    matrix_size=np.shape(Mask)
# extract_CSF.m:17
    
#下一行调试正确后可代替X,Y,Z的定义
#    X,Y,Z=np.mgrid(np.dot((np.arange(0,matrix_size[0])),voxel_size[0]),np.dot((np.arange(0,matrix_size[1])),voxel_size[1]),np.dot((np.arange(0,matrix_size[2])),voxel_size[2]),nargout=3)
# extract_CSF.m:22
    tx=np.array([voxel_size[0]]).T
    ty=np.array([voxel_size[1]]).T
    tz=np.array([voxel_size[2]]).T
    X=np.dot((np.arange(0,matrix_size[0])),tx)
    Y=np.dot((np.arange(0,matrix_size[1])),ty)
    Z=np.dot((np.arange(0,matrix_size[2])),tz)
 #   print(X[Mask>0])
    print(X.shape)
    X_cen=sum(X[Mask > 0]) / sum(ravel(Mask))
# extract_CSF.m:23
    Y_cen=sum(Y(Mask > 0)) / sum(ravel(Mask))
# extract_CSF.m:24
    Z_cen=sum(Z(Mask > 0)) / sum(ravel(Mask))
# extract_CSF.m:25
    # Center region (sphere)
    radius_cen=30
# extract_CSF.m:27
    Mask_cen=sqrt(abs(X - X_cen) ** 2 + abs(Y - Y_cen) ** 2 + abs(Z - Z_cen) ** 2) <= radius_cen
# extract_CSF.m:28
    if flag_erode:
        Mask=SMV(Mask,matrix_size,voxel_size,10) > 0.999
# extract_CSF.m:32
    
    Mask_raw_1=multiply((R2s < thresh_R2s),Mask_cen)
# extract_CSF.m:37
    CC=bwconncomp(Mask_raw_1,6)
# extract_CSF.m:38
    numPixels=cellfun(numel,CC.PixelIdxList)
# extract_CSF.m:39
    numPixels_sorted,idxs=sort(numPixels,2,'descend',nargout=2)
# extract_CSF.m:40
    ROIs_region_cen=zeros(matrix_size)
# extract_CSF.m:41
    for i in arange(1,n_region_cen).reshape(-1):
        numPixels_sorted(i)
        idx=idxs(i)
# extract_CSF.m:44
        ROIs_region_cen[CC.PixelIdxList[idx]]=i
# extract_CSF.m:45
    
    Mask_raw_2=multiply((R2s < thresh_R2s),Mask)
# extract_CSF.m:48
    CC=bwconncomp(Mask_raw_2,6)
# extract_CSF.m:49
    numPixels=cellfun(numel,CC.PixelIdxList)
# extract_CSF.m:50
    numPixels_sorted,idxs=sort(numPixels,2,'descend',nargout=2)
# extract_CSF.m:51
    ROIs_region=zeros(matrix_size)
# extract_CSF.m:52
    for i in arange(1,length(idxs)).reshape(-1):
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
#以下为测试所用
R2s=4
Mask=np.array([[[1,1,4,2,1],[1,5,3,6,2],[1,2,6,4,3],[2,4,3,1,3],[3,3,4,2,5]],
               [[1,2,6,3,4],[3,4,6,2,4],[3,7,6,3,2],[2,4,3,1,3],[3,3,4,2,5]],
               [[3,4,6,2,4],[3,7,6,3,2],[2,6,4,5,3],[2,4,3,1,3],[3,3,4,2,5]],
               [[1,1,4,2,1],[1,5,3,6,2],[1,2,6,4,3],[2,4,3,1,3],[3,3,4,2,5]],
               [[1,1,4,2,1],[1,5,3,6,2],[1,2,6,4,3],[2,4,3,1,3],[3,3,4,2,5]]])
flag_erode=5
thresh_R2s=7
voxel_size=[[3,4,6,2,4],[3,7,6,3,2],[2,6,4,5,3]];
Mask_ROI_CSF = extract_CSF(R2s, Mask, voxel_size, flag_erode, thresh_R2s)
print(Mask_ROI_CSF)