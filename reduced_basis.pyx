#! /usr/bin/python2.7

from __future__ import print_function, division

# Check out signal_model.pyx for explanation
import numpy as np
cimport numpy as np
DTYPE = np.int
ctypedef np.int_t DTYPE_t

from libc.math cimport sin, cos, sqrt

import pyximport; pyximport.install()
import signal_model as sm

def ROMGreedy (np.ndarray params, template, double epsilon):
    """This will generate a reduced basis in the given parameter space.

    Args:

        lambda: A parameter space value. At the moment it is assumed to be of the
            following structure:
                ((param^1_min, param^1_max, param^1_step),
                 (param^2_min, param^2_max, param^2_step),
                 ...
                 (param^n_min, param^n_max, param^n_step))

        template: a function which is used as a basis. This also is used to construct
            the training set.

        epsilon: a measure of the error, when constructing the RB

    Returns:
        
        RB: An array containing the reduced basis parameters. The form is assumed to be
            as follows:
                ((param^1_1, param^1_2, param^1_3, ... , param^1_m),
                 (param^2_1, param^2_2, param^2_3, ... , param^2_m),
                 ...
                 (param^n_1, param^n_2, param^n_3, ... , param^n_m))
    """

    # Initialise some variables
    cdef double projection
    cdef int i = 0, translateIt
    cdef np.ndarray sigma, RB, tmp_l, params_i, params_trial, Grammian
    tmp_l.resize(params.shape[0])
    
    # Initialise the parameter space and calculate the size of the parameter space
    # (total number of points)
    paramSpace = []
    totalNumber = 1
    for i in params:
        paramSpace += [ np.linspace(i[0], i[1], i[2]) ]
        totalNumber = totalNumber * i[2]
    
    # Seed choice (arbitrary)
    RB.resize(len(paramSpace))
    for i in range(RB.shape[0]):
        RB[i] = paramSpace[i][0]

    # Initiate the greedy algorithm
    sigma.resize(1)
    sigma[0] = 1

    # The parameter space is large, so the computation will be expensive
    while sigma[i] > epsilon:
        i += 1
        sigma[i] = 0
        projection = 0
        Grammian = constructGrammian (RB, template)

        # Stupidly traverse the entire parameter space
        # NOTE: We could have MCMC or a simple MC method as well
        for j in range(totalNumber):
            # the params_trial does not have to be zeroed away
            # This is a temporary variable
            translateIt = j

            # Construct the test vector of the parameters
            for k in range(len(paramSpace)):
                n = translateIt % params[k][2]
                params_trial[k] = paramSpace[k][n]
                translateIt = translateIt // params[k][2]

            # Construct the Gram matrix

            # Calculate the projection

            sigmaTrial = np.abs(template(t, params, L, u_p) - projection)
            sigmaTrial *= sigmaTrial

            if sigmaTrial > sigma[i]:
                sigma[i] = sigmaTrial
                params_i = params_trial

        # Add the lambda_i, which was found by maximizing
        RB = np.append(RB, params_i)

    return RB

def constructGrammian (set, template):
    G = np.zeros(set.shape[1], set.shape[1])

    for i in range(G.shape[0]):
        for j in range(G.shape[0]):
            G[i,j] = template(1,2)

    return G
