# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
# .\fgrad.m

    # Discrete Gradient Using Forward Differences
# with the Neuman Boundary Condition
    
    # Created by Youngwook Kee (Oct 21 2015)
# Last modified date: Oct 24 2015
    
    # References:
# [1] Chambolle, An Algorithm for Total Variation Minimization and
# Applications, JMIV 2004
# [2] Pock et al., Global Solutions of Variational Models with Convex
# Regularization, SIIMS 2010
    
def fgrad(chi=None,voxel_size=np.array([1,1,1])):
# .\fgrad.m:16
    
    # chi = double(chi);
    
    Dx=np.concatenate((chi[1:,:,:],chi[-1,:,:])) - chi
# .\fgrad.m:21
    Dy=np.concatenate((chi[:,1:,:],chi[:,-1,:]),axis=1) - chi
# .\fgrad.m:22
    Dz=np.concatenate((chi[:,:,1:],chi[:,:,-1]),axis=2) - chi
# .\fgrad.m:23
    Dx=Dx / voxel_size[0]
# .\fgrad.m:25
    Dy=Dy / voxel_size[1]
# .\fgrad.m:26
    Dz=Dz / voxel_size[2]
# .\fgrad.m:27
    Gx=np.concatenate((Dx,Dy,Dz),axis=3)
# .\fgrad.m:29
    return Gx
    
if __name__ == '__main__':
    pass
    