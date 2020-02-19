# Generated with SMOP  0.41
import numpy as np
import scipy.io as scio

def cgsolve(A,b,tol,maxiter,verbose=None,x0=None):
    
    
    matrix_size = len(b)
    b = b.flatten(1)
    implicit = isinstance(A,type(cgsolve))
    x = np.zeros((len(b),1))
    if x0 == None:
        x = np.zeros((len(b),1))
        r = b  
    else:
        x = x0.reshape(1,-1)
        if implicit:
            r = b-A(np.reshape(x,matrix_size))
            r = r.reshape(1,-1)
        else:
            r = b

    if verbose==None:
        verbose = 1
    
    
    d = r
    delta = r.conj().T*r
    delta = np.array(delta)
    delta0 = (b.conj().T)*b
    numiter = 0
    bestx= x
    bestres = np.sqrt(delta/delta0)
    print('aa',bestres)
    while np.logical_and((numiter < maxiter), (delta> tol**2*delta0) ):
        #q = A*d
        if implicit:
            q = A(np.reshape(d,matrix_size))
            q = q.reshape(1,-1)
        else:
            q = A*d
        alpha = delta/((d.conj().T)*q)
        print(alpha)
        x = x + d.dot(alpha)
        if (numiter+1) % 50 == 0:
            # r = b - Aux*x
            if implicit:
                r = b - np.reshape(A(np.reshape(x,matrix_size)),len(b))
            else:
                r = b - A*x
        else:
            r = r - q.dot(alpha)

        deltaold = delta
        delta = (r.conj().T)*r
        beta = delta / deltaold
        d = r + d.dot(beta)
        numiter = numiter + 1
        item = np.mat(np.sqrt(delta/delta0))
        if (item.any()<bestres.any()):
            bestx = x
            bestres = np.sqrt(delta/delta0)

        if (verbose) and (numiter % verbose) == 0:
            print('cg:Iter = %d' %(numiter))
            print('Best residual = ',bestres)
            print('Current residual = ',np.sqrt(delta/delta0))

    if verbose:
        print('cg:Iterations =',numiter,'best residual = ',bestres)
    #x = reshape(bestx,matrix_size)
    
    x = np.reshape(x,matrix_size)
    res = bestres
    iterr = numiter
    return x,res,iterr