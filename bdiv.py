# Generated with SMOP  0.41
# from smop.libsmop import *
import numpy as np
# .\bdiv.m

    # Discrete Divergence Using Backward Difference
# with the Dirichlet Boundary Condition
    
    # Created by Youngwook Kee (Oct 21 2015)
# Last modified date: Oct 24 2015
    
    # References:
# [1] Chambolle, An Algorithm for Total Variation Minimization and
# Applications, JMIV 2004
    
def bdiv(Gx=None,voxel_size=np.array([1,1,1])):
    
    Gx_x=Gx[:,:,:,0]
    Gx_y=Gx[:,:,:,1]
    Gx_z=Gx[:,:,:,2]
    
    (Mx,My,Mz)=Gx_x.shape[:3]
    
    Dx = np.concatenate((Gx_x[:-1,:,:],np.zeros((1,My,Mz))))        - np.concatenate((np.zeros((1,My,Mz)),Gx_x[:-1,:,:]))
    Dy = np.concatenate((Gx_y[:,:-1,:],np.zeros((Mx,1,Mz))),axis=1) - np.concatenate((np.zeros((Mx,1,Mz)),Gx_y[:,:-1,:]),axis=1)
    Dz = np.concatenate((Gx_z[:,:,:-1],np.zeros((Mx,My,1))),axis=2) - np.concatenate((np.zeros((Mx,My,1)),Gx_z[:,:,:-1]),axis=2)

    Dx=Dx / voxel_size[0]
    Dy=Dy / voxel_size[1]
    Dz=Dz / voxel_size[2]

    div = - (Dx + Dy + Dz)

    return div
    
if __name__ == '__main__':
    pass
    