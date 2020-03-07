import numpy as np

# Discrete Gradient Using Forward Differences
# with the Neuman Boundary Condition
#     
# Created by Youngwook Kee (Oct 21 2015)
# Last modified date: Oct 24 2015
#    
# References:
# [1] Chambolle, An Algorithm for Total Variation Minimization and
# Applications, JMIV 2004
# [2] Pock et al., Global Solutions of Variational Models with Convex
# Regularization, SIIMS 2010
    
def fgrad(chi=None,voxel_size=np.array([1,1,1])):
    # chi = double(chi);
    
    Dx=np.concatenate((chi[1:,:,:],chi[-1:,:,:])) - chi
    Dy=np.concatenate((chi[:,1:,:],chi[:,-1:,:]),axis=1) - chi
    Dz=np.concatenate((chi[:,:,1:],chi[:,:,-1:]),axis=2) - chi

    Dx=Dx / voxel_size[0]
    Dy=Dy / voxel_size[1]
    Dz=Dz / voxel_size[2]
    Gx=np.concatenate((Dx[:,:,:,np.newaxis],Dy[:,:,:,np.newaxis],Dz[:,:,:,np.newaxis]),axis=3)

    return Gx
    
if __name__ == '__main__':
    pass
# chi=np.zeros((256,256,64))
# fgrad(chi)
