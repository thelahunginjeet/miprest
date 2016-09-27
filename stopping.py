#!/usr/bin/env python
# encoding: utf-8

"""
stopping.py

@author: Kevin S. Brown (UCONN), Ameya Akkalkotkar (UCONN)

Created by Kevin Brown on 2016-07-15.
"""
import numpy as np
from numpy import newaxis,dot,log2,log
from numpy.linalg import svd
from scipy.stats import chisqprob

def covmatrix(X):
    '''
    Computes the N x N covariance matrix for an N x p data matrix X.
    '''
    cX = X - X.mean(axis=1)[:,newaxis]
    return dot(cX,cX.T)/(cX.shape[0] - 1)


def corrmatrix(X):
    '''
    Computes the N x N correlation matrix for an N x p data matrix X.
    '''
    sX = (X - X.mean(axis=1)[:,newaxis])/X.std(axis=1)[:,newaxis]
    return dot(sX,sX.T)/(sX.shape[0] - 1)


class StoppingRule(object):
    '''
    Class for computing the number of principal components to retain for a data matrix.
    Covariance and correlation matrices, along with their eigenvalues, are pre-computed
    for the data matrix passed to the class constructor.
    '''
    def __init__(self,data):
        self.data = data
        self.N,self.p = self.data.shape
        self.covMatrix = covmatrix(self.data)
        self.corMatrix = corrmatrix(self.data)
        _,self.covlambda,_ = svd(self.covMatrix, full_matrices=False)
        _,self.corlambda,_ = svd(self.corMatrix, full_matrices=False)


    def bartlett(self, pcrit = 0.05):
        """
        Bartlett's test for the first principal component.  Returns the p-value for
        Bartlett's test statistic.  If the p-value is above the critical value of
        interest, then the data does have even a single principal component.  Therefore,
        Bartlett's test is not a true stopping rule (except insofar as it may say stop
        at zero components).
        """
        b = -(self.N - (2.0*self.p + 11.0)/6.0)*log(self.corlambda).sum()
        dof = self.p*(self.p - 1.0)/2.0
        if chisqprob(b,dof) < pcrit:
            return 1
        return 0

    '''
    def sphericity(self,pcrit = 0.05):
        """
        Sphereicity test of Pimentel (from Bartlett).  Computes the statistic for each
        axis and tests them against pcrit.  The number of axes with the statistic less
        than pcrit is the number to retain.
        """
        s = zeros(self.N)
        for k in xrange(1,len(s)+1):
            indx = k-1
            lbar = self.corlambda[indx+1:]
            s[i] = (self.N - k - (2.0*(self.p - k) + 7.0 + 2.0/(self.p - k))/6.0 + )
    '''

    def kaiser_guttman(self):
        """
        Kaiser Guttman rule selects all components with eigenvalues greater than the average eigenvalue.
        """
        return sum(self.covlambda > self.covlambda.mean())


    def joliffe_kaiser_guttman(self):
        """
        Joliffe's modification to KG which sets the cutoff condition to 0.7*(mean eigenvalue).
        """
        return sum(self.covlambda > 0.7*self.covlambda.mean())


    def broken_stick(self):
        """Broken-stick model: returns a vector of the difference between an eigenvalue and its corresponding
        broken stick counterpart. Drop any signal for which this value is -ve."""
        s0 = self.covlambda
        n = self.N
        l = np.zeros((n,1))
        for j in xrange(0,n):
            l[n-j-1] = (1./n)*np.cumsum(1./(n-j))
        return s0-l.T


    def information_dimension(self):
        """
        Calculates the information dimension (see Cangelosi 2007).
        """
        s_eig = self.covlambda
        pk = self.covlambda/np.sum(self.covlambda)
        h = -1.0*sum(pk*log2(pk))/log2(self.N)
        return self.N**h


    def parallel_analysis(self,pcrit = 0.05):
        """performs the parallel analysis algorithm on the covariance matrix. It creates 1000 random data matrices
        in [0,1] and computes the 95th percentile of random eigenvalues the covariance matrices. Any signal for which
        eig(data)>eig_95pc(random) is retained"""
        s0 = self.covlambda
        s = np.zeros((1000,self.N))
        for i in xrange(0,1000):
            R = np.random.rand(self.N*self.p).reshape((self.N,self.p))
            M0_r = covmatrix(R)
            _,s_temp,_ = np.linalg.svd(M0_r, full_matrices =False)
            s[i,:] = s_temp
        s_crit = np.percentile(s,100*(1-pcrit),axis = 0)
        c = 0
        for i in xrange(len(s0)):
            if s0[i] > s_crit[i]:
                c = c+1
        return c

    '''
    def random_lambda(self,pcrit = 0.05):
        """performs a p-value statistics test on the eigenvalues of 999 permutations of the data matrix"""
        s0 = self.covlambda
        Ng = np.zeros(self.N)
        # compute
        for i in xrange(0,999):
            R = np.random.permutation(self.data.flatten()).reshape((self.data.shape))
            M_r = covmatrix(R)
            _,s_temp,_ = np.linalg.svd(M_r, full_matrices = False)
            Ng += (s_temp > s0)
        #select
        c = 0
        for i in xrange(0,len(Ng)):
            if (Ng[i]+1)/1000.0 < pcrit:
                c = c + 1
        return c
    '''
