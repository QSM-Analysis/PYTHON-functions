# Generated with SMOP  0.41
#from smop.libsmop import *
import numpy as np
import scipy.io as scio
from numpy import sum
# .\parse_QSM_input.m

    # An interal function to parse input arguments for MEDI_L1
#   Created by Tian Liu on 2011.03.16
#   Modified by Tian Liu and Shuai Wang on 2011.03.15
#   Last modified by Tian Liu on 2013.07.24
    
def parse_QSM_input(varargin=None,*args,**kwargs):
    merit=0
# .\parse_QSM_input.m:10
    smv=0
# .\parse_QSM_input.m:11
    lam=1000
# .\parse_QSM_input.m:12
    radius=5
# .\parse_QSM_input.m:13
    data_weighting=1
# .\parse_QSM_input.m:14
    gradient_weighting=1
# .\parse_QSM_input.m:15
    pad=np.array([0])
# .\parse_QSM_input.m:16
    matrix_size0=0
# .\parse_QSM_input.m:17
    Debug_Mode='NoDebug'
# .\parse_QSM_input.m:18
    solver='gaussnewton'
# .\parse_QSM_input.m:19
    # CSF regularization
    lam_CSF=100
# .\parse_QSM_input.m:21
#    filename='RDF.mat'
    filename='medi_siemens_data.mat'
# .\parse_QSM_input.m:22
    if varargin.shape[0] > 0:
        for k in range(varargin.shape[0]):
            if varargin[k].lower()=='filename':
                filename=varargin[k + 1]
# .\parse_QSM_input.m:26
            if varargin[k].lower()=='lambda':
                lam=varargin[k + 1]
# .\parse_QSM_input.m:29
            if varargin[k].lower()=='data_weighting':
                data_weighting=varargin[k + 1]
# .\parse_QSM_input.m:32
            if varargin[k].lower()=='gradient_weighting':
                gradient_weighting=varargin[k + 1]
# .\parse_QSM_input.m:35
            if varargin[k].lower()=='merit':
                merit=1
# .\parse_QSM_input.m:38
            if varargin[k].lower()=='smv':
                smv=1
# .\parse_QSM_input.m:41
                radius=varargin[k + 1]
# .\parse_QSM_input.m:42
            if varargin[k].lower()=='zeropad':
                pad=varargin[k + 1]
# .\parse_QSM_input.m:45
            if varargin[k].lower()=='DEBUG':
                Debug_Mode=varargin[k + 1]
# .\parse_QSM_input.m:48
            if varargin[k].lower()=='lambda_CSF':
                lam_CSF=varargin[k + 1]
# .\parse_QSM_input.m:51
            if varargin[k].lower()=='solver':
                solver=varargin[k + 1]
# .\parse_QSM_input.m:54
    
    #load(filename,'iFreq','RDF','N_std','iMag','Mask','matrix_size','voxel_size','delta_TE','CF','B0_dir')
    data=scio.loadmat(filename)
    iFreq=data['iFreq']
    RDF=data['RDF']
    N_std=data['N_std']
    iMag=data['iMag']
    Mask=data['Mask']
    matrix_size=data['matrix_size']
    voxel_size=data['voxel_size']
    delta_TE=data['delta_TE']
    CF=data['CF']
    B0_dir=np.array(data['B0_dir'].T[0])

    # CSF regularization
#    if ismember('Mask_CSF',who('-file',filename)):
#        Mask_CSF=logical(getfield(load(filename,'Mask_CSF'),'Mask_CSF'))
#    if np.in1d('Mask_CSF',data):
    if data['Mask_CSF'].size != 0
        Mask_CSF=data['Mask_CSF'].astype(bool)
# .\parse_QSM_input.m:62
    else:
        Mask_CSF=np.array([])
# .\parse_QSM_input.m:64

    if not 'delta_TE' in locals().keys():
        delta_TE=input('TE spacing = ')
        with open(filename,'ab') as f:
            scio.savemat(f,{'delta_TE':delta_TE})

    if not 'CF' in locals().keys():
        CF=input('Central Frequency = ')
        with open(filename,'ab') as f:
            scio.savemat(f,{'CF':CF})
    
    if not 'B0_dir' in locals().keys():
        B0_dir=input('B0 direction = ')
        with open(filename,'ab') as f:
            scio.savemat(f,{'B0_dir':B0_dir})
    
    if sum(pad[:]):
        leng=len(pad)
        matrix_size0=matrix_size.copy()
# .\parse_QSM_input.m:83
        matrix_size=matrix_size + pad

# .\parse_QSM_input.m:84
#       iFreq=padarray(iFreq,pad,'post')
        dim=len(iFreq.shape)
        pad1=np.zeros((dim,2))
        for i in range(dim):
            if i>leng:
                pad1[i]=[0,0]
            else:
                pad1[i]=[0,pad[i]]
        iFreq=np.pad(iFreq,pad1,'constant')
# .\parse_QSM_input.m:85
#       RDF=padarray(RDF,pad,'post')
        dim=len(RDF.shape)
        pad1=np.zeros((dim,2))
        for i in range(dim):
            if i>leng:
                pad1[i]=[0,0]
            else:
                pad1[i]=[0,pad[i]]
        iFreq=np.pad(RDF,pad1,'constant')
# .\parse_QSM_input.m:86
#       N_std=padarray(N_std,pad,'post')
        dim=len(N_std.shape)
        pad1=np.zeros((dim,2))
        for i in range(dim):
            if i>leng:
                pad1[i]=[0,0]
            else:
                pad1[i]=[0,pad[i]]
        N_std=np.pad(N_std,pad1,'constant')
# .\parse_QSM_input.m:87
#       iMag=padarray(iMag,pad,'post')
        dim=len(iMag.shape)
        pad1=np.zeros((dim,2))
        for i in range(dim):
            if i>leng:
                pad1[i]=[0,0]
            else:
                pad1[i]=[0,pad[i]]
        iMag=np.pad(iMag,pad1,'constant')
# .\parse_QSM_input.m:88
#       Mask=padarray(Mask,pad,'post')
# .\parse_QSM_input.m:89
        dim=len(Mask.shape)
        pad1=np.zeros((dim,2))
        for i in range(dim):
            if i>leng:
                pad1[i]=[0,0]
            else:
                pad1[i]=[0,pad[i]]
        Mask=np.pad(Mask,pad1,'constant')
        
        if (Mask_CSF.size!=0):
            dim=len(Mask_CSF.shape)
            pad1=np.zeros((dim,2))
            for i in range(dim):
                if i>leng:
                    pad1[i]=[0,0]
                else:
                    pad1[i]=[0,pad[i]]
            Mask_CSF=np.pad(Mask_CSF,pad1,'constant')
# .\parse_QSM_input.m:93
    
    return lam,iFreq,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,radius,data_weighting,gradient_weighting,Debug_Mode,lam_CSF,Mask_CSF,solver
    
if __name__ == '__main__':
    parse_QSM_input()
    
