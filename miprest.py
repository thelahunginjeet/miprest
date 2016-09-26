#!/usr/bin/env python
# encoding: utf-8

"""
miprest.py

@author: Kevin S. Brown (UCONN), Ameya Akkalkotkar (UCONN)

Created by Kevin Brown on 2016-09-19.
"""

from pycar import RAICAR
from rpy2ica import fastica as rfastica
from pyica import fastica
from numpy import ceil,tile,newaxis
import pylab


def decimate_matrix(X,dfactor):
    '''
    Returns a matrix X with X.shape[0] - ceil(X.shape[1]/dfactor) random
    columns removed.  So if dfactor = 4 and X is 10 x 1000, the matrix
    returned will be size 10 x 250.
    '''
    ndrop = X.shape[1] - ceil(X.shape[1]/dfactor)
    dropcol = permutation(X.shape[1])[:ndrop]
    return delete(X,dropcol,1)


class MIPReSt(object):
    '''
    Class for Mixed ICA/PCA via the MIPReSt algorithm.
    '''
    def __init__(self,projDirectory = './scratch',icaOptions=None,avgMethod='selective',icaMethod=fastica):
        '''
        The ProjectDirectory will store intermediate files needed for the
        RAICAR algorithm.  icaOptions should be specified as a dictionary;
        default options are used otherwise.

        INPUT:
        ------
        projDirectory: string, optional
            temporary directory to store ICA calculations during RAICAR

        icaOptions: dict, optional
            see fastica or rfastica for information on ICA keyword arguments; some
            options are not compatible with both methods

        avgMethod: string, optional
            see pycar.RAICAR for details; 'selective' is the original RAICAR
            averaging method.

        icaMethod: callable, optional
            function that performs ICA and returns the mixing matrix, unmixing
            matrix, and the sources.  Used within RAICAR (see the RAICAR package
            for more details).
        '''
        self.projDirectory = projDirectory
        if icaOptions is None:
            self.icaOptions = {}
            self.icaOptions['algorithm'] = 'parallel fp'
            self.icaOptions['nonlinearity'] = 'logcosh'
            self.icaOptions['maxIterations'] = 500
            self.icaOptions['tolerance'] = 1.0e-05
        else:
            self.icaOptions = icaOptions
        self.avgMethod = avgMethod
        self.icaMethod = icaMethod
        return


    def parent_decomposition(self,X,avgMethod='selective'):
        '''
        Performs RAICAR on the parent (undecimated) data set X.  The sources,
        reproducibilities, and mixing matrix are stored in MIPReSt.parent
        '''
        self.parent = {}
        raic = raicar.RAICAR(projDirectory = 'scratch',nSignals=X.shape[0],K=30,avgMethod=self.avgMethod,icaMethod=self.icaMethod,icaOptions=self.icaOptions)
        raic.clean_project()
        raic.runall(X)
        self.parent['R'] = raic.read_reproducibility()
        self.parent['S'],self.parent['A'] = raic.read_raicar_components()
        return


    def decimated_decomposition(self,X,ndec=10,dfactor=2):
        '''
        Performs RAICAR runs on the input data nrep times, decimating by a
        factor of dfactor.  By default, full rank extraction of the matrix
        X (extraction of the same number of signals as the number of rows of
        X) is performed.

        INPUT:
        ------
        X : array, required
            X should be an Nmix x Nsamples matrix of data

        dfactor : integer, optional
            decimation factor; decimated data will have a number of samples
            equal to ceiling(Nsamples/dfactor)

        ndec : integer, optional
            number of random decimations to extract sources from

        OUTPUT:
        ------

        '''
        # structure to hold the results
        self.children = {}.fromkeys(range(ndec))
        # RAICAR object to peform the decompositions
        raic = raicar.RAICAR(projDirectory = 'scratch',nSignals=X.shape[0],K=30,avgMethod=self.avgMethod,icaMethod=self.icaMethod,icaOptions=self.icaOptions)
        # loop over decimations
        for d in range(ndec):
            Xd = decimate_matrix(X,dfactor)
            raic.clean_project()
            raic.runall(Xd)
            self.children[d]['R'] = raic.read_reproducibility()
            # do we need to save these?
            self.children['S'],self.children[d]['A'] = raic.read_raicar_components()
        return


    def compute_R_delta(self):
        '''
        Sorts and packs reproducibilty data into a matrix of size ndecimations
        x nsignals, and computes delta_ij and packs them into a similar shape
        array.  Used for plot_R_delta, in order to determine the dimension of
        the sparse subspace.
        '''
        # check to make sure the parent and child data exists
        if self.parent is None or self.children is None:
            print('You must run both decompositions (parent and decimated) first.')
            return None,None
        # figure out how many decimations and how many signals
        ndec = len(self.children.keys()) + 1    # parent is included!
        nsig = len(self.parent['R'])
        R = zeros((ndec,nsig))
        # transfer R values
        R[0,:] = sort(self.parent['R'])[::-1]
        for i in range(1,ndec):
            R[i,:] = sort(self.children[i]['R'])[::-1]
        # compute delta_ij
        delta[i][j] = abs(R - tile(R[0,:],(ndec,1)))
        return R,deltaij


    def project_sparse(self,X,nsparse):
        '''
        Projects the first nsparse parent RAICAR components out of the data matrix
        X, returning a matrix that consists only of Gaussian components.
        '''
        # sort by reproducibility
        sortedindx = argsort(self.parent['R'])[::-1]
        A = self.parent['A'][:,sortedindx[:nsparse]]
        S = self.parent['S'][sortedindx[:nsparse],:]
        return X - dot(A,S)
