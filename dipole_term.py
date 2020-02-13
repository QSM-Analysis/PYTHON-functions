# Generated with SMOP  0.41
from smop.libsmop import *
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
    

def dipole_term(W=None,D=None,Mask=None,xx=None,*args,**kwargs):
    

    x=np.zeros((D.shape))
# dipole_term.m:17
    #x[(Mask.flatten(1)) == 0]=0
    x[ravel(Mask) == np.zeros(198*212*102)]=np.zeros(D.shape)
# dipole_term.m:18
    x[ravel(Mask) == np.ones(198*212*102)]=np.zeros(D.shape)
# dipole_term.m:19
    Ax=(np.fft.ifftn(multiply(D,np.fft.fftn(multiply(W,(np.fft.ifftn(multiply(D,np.fft.fftn(x)))))))))
# dipole_term.m:21
    y=Ax[ravel(Mask) == np.zeros(198*212*102)]
    print('pass')
# dipole_term.m:22