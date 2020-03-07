import numpy as np
from numpy import max,abs,dot,multiply,sum

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
    
def gradient_mask(iMag=None,Mask=None,grad=None,voxel_size=None,percentage=0.9):
    # modified from matlab code, gradient_weighting_mode is useless, while the percentage is not transfer from GUI
    
    field_noise_level = 0.01 * iMag[:].max()
    wG = abs(grad(iMag * (Mask > 0),voxel_size))
    denominator = sum(Mask[:] == 1)
    numerator = sum(wG[:] > field_noise_level)
    if (numerator / denominator) > percentage:
        while (numerator / denominator) > percentage:
            field_noise_level = field_noise_level * 1.05
            numerator = sum(wG[:] > field_noise_level)

    else:
        while (numerator / denominator) < percentage:
            field_noise_level = field_noise_level * 0.95
            numerator=sum(wG[:] > field_noise_level)
    
    wG = (wG <= field_noise_level)

    return wG