# Generated with SMOP  0.41
from smop.libsmop import *
from bdiv import*
from fgrad import*
from parse_QSM_input import*
from dipole_kernel import*
from SMV import*
import numpy as np
from numpy import multiply,dot,sqrt,abs,exp,mean,pi
from sphere_kernel import*
from dataterm_mask import*
from gradient_mask import*
from cgsolve import*
from store_QSM_results import*
import time
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
    
def MEDI_L1(varargin=None,*args,**kwargs):
    (lambda_,__,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,radius,data_weighting,gradient_weighting,Debug_Mode,lam_CSF,Mask_CSF,solver)=parse_QSM_input(varargin[:])
# .\MEDI_L1.m:43
    ############### weights definition ##############
    print(Mask_CSF.shape)
    cg_max_iter=100
# .\MEDI_L1.m:46
    cg_tol=0.01
# .\MEDI_L1.m:47
    max_iter=10
# .\MEDI_L1.m:48
    tol_norm_ratio=0.1
# .\MEDI_L1.m:49
    data_weighting_mode=data_weighting
# .\MEDI_L1.m:50
    gradient_weighting_mode=gradient_weighting
# .\MEDI_L1.m:51
    grad=fgrad
# .\MEDI_L1.m:52
    div=bdiv
# .\MEDI_L1.m:53
    # grad = @cgrad;
# div = @cdiv;
    
    N_std=multiply(N_std,Mask)
# .\MEDI_L1.m:57
    tempn=N_std.astype(float)
# .\MEDI_L1.m:58
    D=dipole_kernel(np.array([matrix_size,voxel_size,B0_dir]))
# .\MEDI_L1.m:59
    if (smv):
        #     S = SMV_kernel(matrix_size, voxel_size,radius);
        SphereK=sphere_kernel(matrix_size,voxel_size,radius).astype(float)
# .\MEDI_L1.m:63
        Mask=SMV(Mask,[SphereK]) > 0.999
# .\MEDI_L1.m:64
        D=multiply((1 - SphereK),D)
# .\MEDI_L1.m:65
        RDF=RDF - SMV(RDF,[SphereK])
# .\MEDI_L1.m:66
        RDF=multiply(RDF,[Mask])
# .\MEDI_L1.m:67
        tempn=sqrt(SMV(tempn ** 2,[SphereK]) + tempn ** 2)
# .\MEDI_L1.m:68
    
    m=dataterm_mask(data_weighting_mode,tempn,Mask)
# .\MEDI_L1.m:71
    b0=multiply(m,exp(dot(1j,RDF)))
# .\MEDI_L1.m:72
    wG=gradient_mask(gradient_weighting_mode,iMag,Mask,grad,voxel_size)
# .\MEDI_L1.m:73
    # CSF regularization
    flag_CSF=(Mask_CSF.size!=0)
# .\MEDI_L1.m:77
    if flag_CSF:
        print('CSF regularization used\n')
    
    oldN_std=N_std.copy()
# .\MEDI_L1.m:81
    print('Using ',solver,'\n')

# .\MEDI_L1.m:85
    def gaussnewton(lambda_,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,radius,Debug_Mode,lam_CSF,Mask_CSF,
    	cg_tol,cg_max_iter,max_iter,tol_norm_ratio,data_weighting_mode,gradient_weighting_mode,grad,div,tempn,D,m,b0,wG,flag_CSF,*args,**kwargs):

        if flag_CSF:
            LT_reg=lambda x=None: multiply(Mask_CSF,(x - mean(x[Mask_CSF])))
    # .\MEDI_L1.m:91
           
        iter=-1
    # .\MEDI_L1.m:94
        x=np.zeros(matrix_size.astype(int)[0])
    # .\MEDI_L1.m:95
        
        # if (logical_not(isempty(findstr(upper(Debug_Mode),'SAVEITER')))):
        #     store_CG_results(multiply(dot(x / (dot(dot(dot(2,pi),delta_TE),CF)),1000000.0),Mask))
        
        res_norm_ratio=np.inf
    # .\MEDI_L1.m:99
        cost_data_history=np.zeros((1,max_iter))
    # .\MEDI_L1.m:100
        cost_reg_history=np.zeros((1,max_iter))
    # .\MEDI_L1.m:101
        e=1e-06
    # .\MEDI_L1.m:103
        
        badpoint=np.zeros(matrix_size.astype(int)[0])
    # .\MEDI_L1.m:104
        Dconv=lambda dx=None: np.real(np.fft.ifftn(multiply(D,np.fft.fftn(dx))))
    # .\MEDI_L1.m:105
        while (res_norm_ratio > tol_norm_ratio) and (iter < max_iter):

            tic()
            iter=iter + 1
    # .\MEDI_L1.m:108
            Vr=1.0 / sqrt(abs(multiply(wG,grad(np.real(x),voxel_size))) ** 2 + e)
    # .\MEDI_L1.m:109
            w=multiply(m,exp(dot(1j,np.fft.ifftn(multiply(D,np.fft.fftn(x))))))
    # .\MEDI_L1.m:110
            #reg=lambda dx=None: div(multiply(wG,(multiply(Vr,(multiply(wG,grad(np.real(dx),voxel_size)))))),voxel_size)
            reg = lambda dx=None: div(wG*(Vr*(wG*grad(np.real(dx),voxel_size))),voxel_size)
            reg1=lambda dx=None: div(multiply(wG,(multiply(Vr,(multiply(wG,grad(np.real(dx),voxel_size)))))),voxel_size)
    # .\MEDI_L1.m:111
            if flag_CSF:
                reg_CSF=lambda dx=None: multiply(lam_CSF,LT_reg(LT_reg(np.real(dx))))
    # .\MEDI_L1.m:113
                reg=lambda dx=None: reg1(dx) + reg_CSF(dx)
    # .\MEDI_L1.m:114
            fidelity=lambda dx=None: Dconv(multiply(multiply(w.conj(),w),Dconv(dx)))
    # .\MEDI_L1.m:116
            A=lambda dx=None: reg(dx) + dot(dot(2,float(lambda_)),fidelity(dx))
    # .\MEDI_L1.m:118
            b=reg(x) + dot(dot(2,float(lambda_)),Dconv(np.real(multiply(multiply(w.conj(),-1j),(w-b0)))))
    # .\MEDI_L1.m:119
            (dx,res,iterr)=np.real(cgsolve(A,-b[0],cg_tol,cg_max_iter,0))
    # .\MEDI_L1.m:123
            res_norm_ratio=np.linalg.norm(ravel(dx)) / np.linalg.norm(ravel(x))
    # .\MEDI_L1.m:124
            x=x + dx
    # .\MEDI_L1.m:125
            wres=multiply(m,exp(dot(1j,(np.real(np.fft.ifftn(multiply(D,np.fft.fftn(x)))))))) - b0
    # .\MEDI_L1.m:127
            cost_data_history[0][iter]=np.linalg.norm(ravel(wres))
    # .\MEDI_L1.m:129
            cost=abs(multiply(wG,grad(x)))
    # .\MEDI_L1.m:130
            cost_reg_history[0][iter]=sum(ravel(cost))
    # .\MEDI_L1.m:131
            if merit:
                wres=wres - mean(wres[ravel(Mask) == 1])
    # .\MEDI_L1.m:135
                a=wres[ravel(Mask) == 1]
    # .\MEDI_L1.m:136
                factor=dot(np.std(abs(a),ddof=1),6)
    # .\MEDI_L1.m:137
                wres=abs(wres) / factor
    # .\MEDI_L1.m:138
                wres[wres < 1]=1
    # .\MEDI_L1.m:139
                badpoint[wres > 1]=1
    # .\MEDI_L1.m:140
                N_std[Mask == 1]=multiply(N_std[Mask == 1],wres[Mask == 1] ** 2)
    # .\MEDI_L1.m:141
                tempn=N_std.astype(float)
    # .\MEDI_L1.m:142
                if (smv):
                    tempn=sqrt(SMV(tempn ** 2,SphereK) + tempn ** 2)
    # .\MEDI_L1.m:144
                m=dataterm_mask(data_weighting_mode,tempn,Mask)
    # .\MEDI_L1.m:146
                b0=multiply(m,exp(dot(1j,RDF)))
    # .\MEDI_L1.m:147
            print('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost:%8.4f.\n' %(iter,res_norm_ratio,cost_data_history[0][iter],cost_reg_history[0][iter]))
            toc()  
        
        #convert x to ppm
        x=multiply(dot(x / (dot(dot(dot(2,pi),float(delta_TE)),CF)),1000000.0),Mask)
    # .\MEDI_L1.m:159
        
        if flag_CSF:
            x=x - mean(x[Mask_CSF])
    # .\MEDI_L1.m:163
               
        if (matrix_size0):
            x=x[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
    # .\MEDI_L1.m:167
            iMag=iMag[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
    # .\MEDI_L1.m:168
            RDF=RDF[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
    # .\MEDI_L1.m:169
            Mask=Mask[:matrix_size0[0], :matrix_size0[1], :matrix_size0[2]]
    # .\MEDI_L1.m:170
            matrix_size=matrix_size0.copy()
    # .\MEDI_L1.m:171
        

    #   resultsfile=store_QSM_results(x,iMag,RDF,Mask,'Norm','L1','Method','MEDIN','Lambda',lambda_,'SMV',smv,'Radius',radius,'IRLS',merit,'voxel_size',voxel_size,'matrix_size',matrix_size,'Data_weighting_mode',data_weighting_mode,'Gradient_weighting_mode',gradient_weighting_mode,'L1_tol_ratio',tol_norm_ratio,'Niter',iter,'CG_tol',cg_tol,'CG_max_iter',cg_max_iter,'B0_dir',B0_dir)
# .\MEDI_L1.m:174
        
        return x,cost_reg_history,cost_data_history
    if 'gaussnewton' == solver:
        x,cost_reg_history,cost_data_history=gaussnewton(lambda_,RDF,N_std,iMag,Mask,matrix_size,matrix_size0,voxel_size,delta_TE,CF,B0_dir,merit,smv,radius,Debug_Mode,lam_CSF,Mask_CSF,
    	cg_tol,cg_max_iter,max_iter,tol_norm_ratio,data_weighting_mode,gradient_weighting_mode,grad,div,tempn,D,m,b0,wG,flag_CSF)
    resultsfile=''
    return x,cost_reg_history, cost_data_history, resultsfile
if __name__ == '__main__':
    pass
