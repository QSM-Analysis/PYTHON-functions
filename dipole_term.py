# Generated with SMOP  0.41

from sklearn.externals import joblib
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
    W=joblib.load('W.pkl')
    D=joblib.load('D.pkl')
    Mask=joblib.load('Mask.pkl')

    x=np.zeros((D.shape))
# dipole_term.m:17
    #x[(Mask.flatten(1)) == 0]=0
    x[Mask == 0]=xx
# dipole_term.m:18
    x[Mask == 1]=0
# dipole_term.m:19
    Ax=np.real(np.fft.ifftn(D*np.fft.fftn(W*np.real(np.fft.ifftn(D*np.fft.fftn(x))))))
# dipole_term.m:21
    y=Ax[Mask == 0]

    return y
    
# dipole_term.m:22