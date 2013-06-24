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

def ROMGreedy (np.ndarray lambda, h, double epsilon):
    """This will generate a reduced basis in the given parameter space.

    Args:

        lambda: A parameter space value. At the moment it is assumed to be of the
            following structure:
                ((param^1_min, param^1_max, param^1_step),
                 (param^2_min, param^2_max, param^2_step),
                 ...
                 (param^n_min, param^n_max, param^n_step))

        h: a function which is used as a basis. This also is used to construct the
            training set.

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
    cdef int i = 0
    cdef np.ndarray sigma, RB, tmp_l
    tmp_l.resize(lambda.shape[0])
    
    # Initialise the parameter space
    paramSpace = []
    for i in lambda:
        paramSpace += [ np.linspace(i[0], i[1], i[2]) ]

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

        # Find the 
        # find the max of sigma = || h - p(i-1)h||^2
        # Save the lambda value when sigma is max

        # RB = np.append(RB, lambda_i)

    return RB
