# from sklearn.externals import joblib
import numpy as np

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
    

def dipole_term(W,D,Mask,xx):
    x=np.zeros((D.shape))
    x[Mask == 0]=xx
    x[Mask == 1]=0
    
    Ax=np.real(np.fft.ifftn(D*np.fft.fftn(W*np.real(np.fft.ifftn(D*np.fft.fftn(x))))))
    y=Ax[Mask == 0]

    return y
    
# dipole_term.m:22
