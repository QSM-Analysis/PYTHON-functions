# Generated with SMOP  0.41
from smop.libsmop import *
import PDF
import numpy as np
# dipole_term.m

    # Normal Equation of the Forward Calculation used in PDF
#   y = dipole_term(W,D,Mask,xx)
# 
#   output
#   y - a background field 
# 
#   input
#   W - noise covariance matrix
#   D - dipole kernel
#   xx - the background dipoles
    
    #   Created by Tian Liu in 2009
#   Last modified by Tian Liu on 2013.07.24
    

def dipole_term(xx,*args,**kwargs):
    
    
    x=np.zeros((PDF.D_temp.shape))
# dipole_term.m:17
    #x[(Mask.flatten(1)) == 0]=0
    x[PDF.Mask_temp == 0]=xx
# dipole_term.m:18
    x[PDF.Mask_temp == 1]=0
# dipole_term.m:19
    Ax=(np.fft.ifftn((PDF.D_temp*np.fft.fftn((PDF.W_temp*(np.fft.ifftn((PDF.D_temp*np.fft.fftn(x)))))))))
# dipole_term.m:21
    y=Ax[PDF.Mask_temp == 0]
    print('pass')
    return y
    
# dipole_term.m:22