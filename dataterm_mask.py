# Generated with SMOP  0.41
from smop.libsmop import *
import numpy as np
from np import multiply,mean
# .\dataterm_mask.m

    # Generate the data weighting
#   w=dataterm_mask(dataterm_weighting_mode, N_std, Mask,cutoff)
# 
#   output
#   w - data weighting
    
    #   input
#   dataterm_weighting_mode - 0, uniform weighting; 1, SNR weighting
#   N_std - noise standard deviation on the field map
#   Mask is a binary 3D matrix denoting the Region Of Interest
    
    #   Created by Ildar Khalidov in 20010
#   Last modified by Tian Liu on 2013.07.24
    
    
@function
def dataterm_mask(dataterm_weighting_mode=None,N_std=None,Mask=None,*args,**kwargs):

    if 0 == dataterm_weighting_mode:
        w=1
# .\dataterm_mask.m:21
    else:
        if 1 == dataterm_weighting_mode:
            w=Mask / N_std
# .\dataterm_mask.m:23
            w[np.isnan(w)]=0
# .\dataterm_mask.m:24
            w[np.isinf(w)]=0
# .\dataterm_mask.m:25
            w=multiply(w,(Mask > 0))
# .\dataterm_mask.m:26
            w=w / mean(w[Mask > 0])
# .\dataterm_mask.m:27
    