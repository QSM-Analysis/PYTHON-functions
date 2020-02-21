# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
from numpy import*
# arlo.m

    
@function
def arlo(te,y,*args,**kwargs):
    varargin = arlo.varargin
    nargin = arlo.nargin

    # r2 = arlo(te,y)
#    
# output
#   r2 r2-map (in Hz)
#      empty when only one echo is provided
    
    # input
#   te  array containing te values (in s)
#   y   a multi-echo data set of arbitrary dimension
#       echo should be the last dimension
#       
# If you use this function please cite
# 
# Pei M, Nguyen TD, Thimmappa ND, Salustri C, Dong F, Cooper MA, Li J, 
# Prince MR, Wang Y. Algorithm for fast monoexponential fitting based 
# on Auto-Regression on Linear Operations (ARLO) of data. 
# Magn Reson Med. 2015 Feb;73(2):843-50. doi: 10.1002/mrm.25137. 
# Epub 2014 Mar 24. PubMed PMID: 24664497; 
# PubMed Central PMCID:PMC4175304.
    
    nte=len(te)
# arlo.m:22
    if nte < 2:
        r2=[]
# arlo.m:24
        return r2
    
    sz=np.shape(y)
# arlo.m:28
    edx=np.size(sz)
# arlo.m:29
    if sz[edx-1] !=nte:
        np.error(np.concat(('Last dimension of y has size ',np.num2str(sz(edx)),', expected ',np.num2str(nte))))
    
    yy=np.zeros((np.arange(0,edx-2)))
# arlo.m:35
    yx=np.zeros((np.arange(0,edx-2)))
# arlo.m:36
    beta_yx=np.zeros((np.arange(0,edx-2)))
# arlo.m:37
    beta_xx=np.zeros((np.arange(0,edx-2)))
# arlo.m:38
    s1=[]
# arlo.m:39
    d1=[]
# arlo.m:39
    crd=tile(':',(0,edx-1))
# arlo.m:40
    crd0=copy(crd)
# arlo.m:41
    crd1=copy(crd)
# arlo.m:41
    crd2=copy(crd)
# arlo.m:41
    for j in arange(0,nte-3).reshape(-1):
        alpha=dot((te[j + 2] - te[j]),(te[j + 2] - te[j])) / 2 / (te[j + 1] - te[j])
# arlo.m:43
        tmp=(dot(dot(2,te[j + 2]),te[j + 2]) - dot(te[j],te[j + 2]) - dot(te[j],te[j]) + dot(dot(3,te[j]),te[j + 1]) - dot(dot(3,te[j + 1]),te[j + 2])) / 6
# arlo.m:44
        beta=tmp / (te[j + 2] - te[j + 1])
# arlo.m:45
        gamma=tmp / (te[j + 1] - te[j])
# arlo.m:46
        crd0[edx-1]=j
# arlo.m:47
        crd1[edx-1]=j + 1
# arlo.m:47
        crd2[edx-1]=j + 2
# arlo.m:47
        #     [te(j+2)-te(j)-alpha+gamma alpha-beta-gamma beta]/((te(2)-te(1))/3)
        y1=dot(y(crd0[arange()]),(te(j + 2) - te(j) - alpha + gamma)) + dot(y(crd1[arange()]),(alpha - beta - gamma)) + dot(y(crd2[arange()]),beta)
# arlo.m:49
        x1=y(crd0[arange()]) - y(crd2[arange()])
# arlo.m:50
        yy=yy + multiply(y1,y1)
# arlo.m:51
        yx=yx + multiply(y1,x1)
# arlo.m:52
        beta_yx=beta_yx + multiply(dot(beta,y1),x1)
# arlo.m:53
        beta_xx=beta_xx + multiply(dot(beta,x1),x1)
# arlo.m:54
    
    r2=(yx + beta_xx) / (beta_yx + yy)
# arlo.m:57
    r2[isnan(r2)]=0
# arlo.m:58
    r2[isinf(r2)]=0
# arlo.m:59
te=[2,6,2,5,4];
y=[[3,4,6,2,4],[3,7,6,3,2],[2,6,4,5,3]];
r2=arlo(te,y)
print(r2)