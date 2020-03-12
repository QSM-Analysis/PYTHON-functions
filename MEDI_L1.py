# -*- coding: utf-8 -*-
# from smop.libsmop import *
# from bdiv import*
# from fgrad import*
# from parse_QSM_input import*
# from dipole_kernel import*
# from SMV import*
# from sphere_kernel import*
# from dataterm_mask import*
# from gradient_mask import*
# from cgsolve import*
# from store_QSM_results import*

import time
import numpy as np
from numpy.fft import fftn, ifftn
from numpy import multiply,dot,sqrt,abs,exp,mean,pi, ravel
from PyQSM.fgrad import fgrad
from PyQSM.bdiv import bdiv
from PyQSM.dipole_kernel import dipole_kernel
from PyQSM.sphere_kernel import sphere_kernel
from PyQSM.SMV import SMV
from PyQSM.dataterm_mask import dataterm_mask
from PyQSM.cgsolve import cgsolve
from PyQSM.parse_QSM_input import parse_QSM_input
from PyQSM.gradient_mask import gradient_mask
from PyQSM.extract_CSF import isempty
def tic():
    globals()['tt'] = time.clock()
def toc():
    print ('\nElapsed time: %.8f seconds\n' % (time.clock()-globals()['tt']))

# .\MEDI_L1.m

    # Morphology Enabled Dipole Inversion (MEDI)
#   [x, cost_reg_history, cost_data_history] = MEDI_L1(varargin)
    
    #   output
#   x - the susceptibility distribution 
#   cost_reg_history - the cost of the regularization term
#   cost_data_history - the cost of the data fidelity term
#   
#   input
#   RDF.mat has to be in current folder.  
#   MEDI_L1('lambda',lam,...) - lam specifies the regularization parameter
#                               lam is in front of the data fidelity term
    
    #   ----optional----   
#   MEDI_L1('smv', radius,...) - specify the radius for the spherical mean
#                                value operator using differential form
#   MEDI_L1('merit',...) - turn on model error reduction through iterative
#                          tuning
#   MEDI_L1('zeropad',padsize,...) - zero pad the matrix by padsize
#   MEDI_L1('lambda_CSF',lam_CSF,...) - automatic zero reference (MEDI+0)
#                                       also require Mask_CSF in RDF.mat
    
    #   When using the code, please cite 
#   Z. Liu et al. MRM 2017;DOI:10.1002/mrm.26946
#   T. Liu et al. MRM 2013;69(2):467-76
#   J. Liu et al. Neuroimage 2012;59(3):2560-8.
#   T. Liu et al. MRM 2011;66(3):777-83
#   de Rochefort et al. MRM 2010;63(1):194-206
    
    #   Adapted from Ildar Khalidov
#   Modified by Tian Liu on 2011.02.01
#   Modified by Tian Liu and Shuai Wang on 2011.03.15
#   Modified by Tian Liu and Shuai Wang on 2011.03.28 add voxel_size in grad and div
#   Last modified by Tian Liu on 2013.07.24
#   Last modified by Tian Liu on 2014.09.22
#   Last modified by Tian Liu on 2014.12.15
#   Last modified by Zhe Liu on 2017.11.06
    
def MEDI_L1(RDF,N_std,iMag,Mask,Mask_CSF, voxel_size,B0_dir,CF, delta_TE,   # data need input
            _lambda=1000, edge_percentage=0.9, smv_radius=5,                # data change data looking
            max_iter=10,tol_norm_ratio=0.1, cg_max_iter=100, cg_tol=0.01):  # data for iteration and cal time
    #(_lambda,__,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,radius,data_weighting,gradient_weighting,Debug_Mode,lam_CSF,Mask_CSF,solver)=parse_QSM_input(varargin[:])
    # matrix_size,voxel_size = matrix_size[0], voxel_size.ravel()

    # parameter no need to set as input
    matrix_size = np.array(Mask.shape).astype(float)
    if smv_radius ==0: 
        smv = 0
    else:
        smv = 1
    
    # parameter not need to adjust from input, may be use in future
    # pad=np.array([0])
    matrix_size0=0 # zero pad size, assume no pad, pad implement not prepared, consulting original matlab code
    merit = 0
    Debug_Mode='NoDebug'
    solver='gaussnewton'
    lam_CSF=100
    
    ############### weights definition ##############
#     cg_max_iter=100
#     cg_tol=0.01
#     max_iter=1 # 10
#     tol_norm_ratio=0.1
    data_weighting_mode = 1
    gradient_weighting_mode = 1 # useless parameter
    grad=fgrad
    div=bdiv
    
    N_std= N_std * Mask
    tempn=N_std.astype(float)
    D=dipole_kernel(matrix_size,voxel_size,B0_dir)

    if (smv):
        SphereK = sphere_kernel(matrix_size,voxel_size,smv_radius)
        Mask = SMV(Mask.astype(float),[SphereK]) > 0.999
        D = (1 - SphereK) * D
        RDF = RDF - SMV(RDF,[SphereK])
        RDF = RDF * Mask
        tempn = sqrt(SMV(tempn ** 2,[SphereK])+ tempn ** 2).astype(float)
    
    m = dataterm_mask(data_weighting_mode,tempn,Mask).astype(float)
    b0 = m * exp(1j * RDF)
    wG = gradient_mask(iMag,Mask,grad,voxel_size,edge_percentage)

    # CSF regularization
    flag_CSF= not isempty(Mask_CSF) # (Mask_CSF.size!=0)
    if flag_CSF:
        print('CSF regularization used')
    else:
        print('CSF not used')
    oldN_std = N_std.copy()
    print('Using ',solver)

    if 'gaussnewton' == solver:
        x,cost_reg_history,cost_data_history = gaussnewton(_lambda,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,smv_radius,Debug_Mode,lam_CSF,Mask_CSF,
        cg_tol,cg_max_iter,max_iter,tol_norm_ratio,data_weighting_mode,gradient_weighting_mode,grad,div,tempn,D,m,b0,wG,flag_CSF,SphereK)
    
    return x,cost_reg_history, cost_data_history

def gaussnewton(_lambda,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,smv_radius,Debug_Mode,lam_CSF,Mask_CSF,
    cg_tol,cg_max_iter,max_iter,tol_norm_ratio,data_weighting_mode,gradient_weighting_mode,grad,div,tempn,D,m,b0,wG,flag_CSF,SphereK):

    if flag_CSF:
        LT_reg=lambda x=None: Mask_CSF * (x - x[Mask_CSF].mean())
       
    _iter = 0
    x = np.zeros(matrix_size.astype(int))
    
    res_norm_ratio=np.inf
    cost_data_history = np.zeros(max_iter)
    cost_reg_history  = np.zeros(max_iter)

    e = 1e-06 # a very small number to avoid /0
    badpoint = np.zeros(matrix_size.astype(int))
    Dconv = lambda dx=None: np.real(ifftn(D * fftn(dx)) )
    while (res_norm_ratio > tol_norm_ratio) and (_iter < max_iter):
        tic()
        
        Vr = 1.0 / sqrt(abs(wG * grad(np.real(x),voxel_size)) ** 2 + e)
        w = m * exp(1j * ifftn(D * fftn(x)))
        reg = lambda dx=None: div(wG * (Vr * (wG * grad(np.real(dx),voxel_size))),voxel_size)
        # reg1= lambda dx=None: div(multiply(wG,(multiply(Vr,(multiply(wG,grad(np.real(dx),voxel_size)))))),voxel_size)
        if flag_CSF:
            reg_CSF = lambda dx=None: lam_CSF * LT_reg(LT_reg(np.real(dx)))
            reg1 = lambda dx=None: reg(dx) + reg_CSF(dx)
        else:
            reg1 = reg   
        fidelity=lambda dx=None: Dconv(w.conj() * w * Dconv(dx))

        A = lambda dx=None: reg1(dx) + 2 * float(_lambda) * fidelity(dx)
        b = reg1(x) + 2 * float(_lambda) * Dconv(np.real(w.conj() * -1j * (w-b0)))
        # b = reg(x) + reg_CSF(x)
        (dx,res,iterr) = np.real(cgsolve(A,-b,cg_tol,cg_max_iter,0))
        res_norm_ratio = np.linalg.norm(dx.ravel()) / np.linalg.norm(x.ravel())
        x = x + dx
        
        wres = m * exp(1j * (np.real(ifftn(D * fftn(x))))) - b0

        cost_data_history[_iter] = np.linalg.norm(wres.ravel())
        cost = abs(wG * grad(x))
        cost_reg_history[_iter] = sum(cost.ravel())

        if merit:
            wres = wres - mean(wres[Mask.ravel() == 1])
            a = wres[Mask.ravel() == 1]
            factor = np.std(abs(a),ddof=1) * 6
            wres = abs(wres) / factor
            wres[wres < 1]=1
            badpoint[wres > 1]=1
            N_std[Mask == 1] = N_std[Mask == 1] * (wres[Mask == 1] ** 2)
            tempn=N_std.astype(float)
            if (smv):
                tempn = sqrt(SMV(tempn ** 2,[SphereK]) + tempn ** 2)
            m = dataterm_mask(data_weighting_mode, tempn, Mask)
            b0 = m * exp(1j * RDF)

        print('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost:%8.4f.' %(_iter,res_norm_ratio,cost_data_history[_iter],cost_reg_history[_iter]))
        
        _iter = _iter + 1
        toc()  
    
    #convert x to ppm
    x = x / (2 *pi * float(delta_TE) * CF) * 1e6 * Mask

    # Zero reference using CSF       
    if flag_CSF:
        x=x - mean(x[Mask_CSF])
           
    if (matrix_size0):
        x   =     x[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
        iMag=  iMag[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
        RDF =   RDF[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
        Mask = Mask[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
        matrix_size=matrix_size0.copy()

#   resultsfile=store_QSM_results(x,iMag,RDF,Mask,'Norm','L1','Method','MEDIN','Lambda',_lambda,'SMV',smv,'Radius',smv_radius,'IRLS',merit,'voxel_size',voxel_size,'matrix_size',matrix_size,'Data_weighting_mode',data_weighting_mode,'Gradient_weighting_mode',gradient_weighting_mode,'L1_tol_ratio',tol_norm_ratio,'Niter',iter,'CG_tol',cg_tol,'CG_max_iter',cg_max_iter,'B0_dir',B0_dir)
    
    return x,cost_reg_history,cost_data_history

if __name__ == '__main__':
    pass


