# Generated with SMOP  0.41
import numpy as np


def cgsolve(A,b,tol,maxiter,verbose=None,x0=None):
    
    matrix_size = b.shape
    print(matrix_size)
    b = b.flatten()
    
#     implicit = (type(A[0])==type(cgsolve))
#     if implicit:
#         fun, W,D,Mask = A
    implicit = (type(A)==type(cgsolve))
        
    if x0 == None:
        x = np.zeros((len(b),1))
        r = b  
    else:
        x = x0.flatten()
        if implicit:
            #r = b-fun(W, D, Mask, np.reshape(x,matrix_size))
            r = b-A(np.reshape(x,matrix_size))
            r = r.reshape(1,-1)
        else:
            r = b

    if verbose==None:
        verbose = 1
    
    d = r
    delta = np.dot(r.conj().T,r)
    delta0 = np.dot((b.conj().T),b)
    numiter = 0
    bestx= x
    bestres = np.sqrt(delta/delta0)
    
    while (numiter < maxiter) and (delta> tol**2 * delta0):
        
        #q = A*d
        if implicit:
            #q = fun(W, D, Mask,np.reshape(d,matrix_size))
            q = A(np.reshape(d,matrix_size)).flatten()
            #q = q.flatten()
            
        else:
            q = A*d
        alpha = delta/np.dot((d.conj().T),q)
        x = x+d.reshape([d.size,1])*alpha
        
        if (numiter+1) % 50 == 0:
            # r = b - Aux*x
            if implicit:
                #r = b - np.reshape(fun(W, D, Mask,np.reshape(x,matrix_size)),len(b))
                r = b - A(np.reshape(x,matrix_size)).flatten()
            else:
                r = b - np.dot(A,x)
        else:
            r = r - alpha*q

        deltaold = delta
        delta = np.dot((r.conj().T),r)
        beta = delta / deltaold
        d = r + beta*d
        numiter = numiter + 1
       
        if (np.sqrt(delta/delta0) < bestres):
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
