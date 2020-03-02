# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
from numpy import max,abs,dot,multiply,sum
# .\gradient_mask.m

    # Generate the gradient weighting in MEDI
#   w=gradient_mask(gradient_weighting_mode, iMag, Mask, grad, voxel_size, percentage)
# 
#   output
#   w - gradient weighting
    
    #   input
#   gradient_weighting_mode - 1, binary weighting; other, reserved for
#                             grayscale weighting
#   iMag - the anatomical image
#   Mask - a binary 3D matrix denoting the Region Of Interest
#   grad - function that takes gradient
#   voxel_size - the size of a voxel
#   percentage(optional) - percentage of voxels considered to be edges.
    
    #   Created by Ildar Khalidov in 20010
#   Modified by Tian Liu and Shuai Wang on 2011.03.28 add voxel_size in grad
#   Modified by Tian Liu on 2011.03.31
#   Last modified by Tian Liu on 2013.07.24
    
def gradient_mask(gradient_weighting_mode=None,iMag=None,Mask=None,grad=None,voxel_size=None,percentage=0.9,*args,**kwargs):

# .\gradient_mask.m:24
    
    field_noise_level=dot(0.01,max(iMag[:]))
# .\gradient_mask.m:29
    wG=abs(grad(multiply(iMag,(Mask > 0)),voxel_size))
# .\gradient_mask.m:30
    denominator=sum(Mask[:] == 1)
# .\gradient_mask.m:31
    numerator=sum(wG[:] > field_noise_level)
# .\gradient_mask.m:32
    if (numerator / denominator) > percentage:
        while (numerator / denominator) > percentage:

            field_noise_level=dot(field_noise_level,1.05)
# .\gradient_mask.m:35
            numerator=sum(wG[:] > field_noise_level)
# .\gradient_mask.m:36

    else:
        while (numerator / denominator) < percentage:

            field_noise_level=dot(field_noise_level,0.95)
# .\gradient_mask.m:40
            numerator=sum(wG[:] > field_noise_level)
# .\gradient_mask.m:41

    
    wG=(wG <= field_noise_level)
# .\gradient_mask.m:45
    return wG