#!/usr/bin/env python
# encoding: utf-8

"""
miprest.py

@author: Kevin S. Brown (UCONN), Ameya Akkalkotkar (UCONN)

Created by Kevin Brown on 2016-09-19.
"""

from stopping import covmatrix
from numpy.linalg import svd
from numpy import dot,newaxis

def pca(X,k):
    '''
    PCA decomposition of matrix X.  X is assumed to be N x p, where p is the
    number of samples (backwards from many PCA implementations).  If you want
    the p x N version, just transpose what comes out of this function.

    k is the number of components to retain  (probably determined by some PCA stopping rule).

    Returns the matrix of eigenvectors of X (the "mixing matrix") and the "signals"
    (projection of the data onto the first k components).
    '''
    # row center the data matrix
    cX = X - X.mean(axis=1)[:,newaxis]
    C = covmatrix(cX)
    # singular value decomp
    _,s,W = svd(C)
    # select first k columns
    W = W[:,:k]
    # compute signal matrix
    S = dot(W.T,X)
    # need to do something about the units
    return W,S
