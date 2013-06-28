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

def ROMGreedy (schedule, params, pulsars, matrix, epsilon):
    """This will generate a reduced basis in the given parameter space.

    Args:

        params: A parameter space value. At the moment it is assumed to be of the
            following structure:
                ((param^1_min, param^1_max, param^1_step),
                 (param^2_min, param^2_max, param^2_step),
                 ...
                 (param^n_min, param^n_max, param^n_step))

        schedule: An array which contains all the ids of the pulsars and when they
        were/going to be measured.

        epsilon: a measure of the error, when constructing the RB

    Returns:
        
        RB: An array containing the reduced basis parameters. The form is assumed to be
            as follows:
                ((param^1_1, param^1_2, param^1_3, ... , param^1_m),
                 (param^2_1, param^2_2, param^2_3, ... , param^2_m),
                 ...
                 (param^n_1, param^n_2, param^n_3, ... , param^n_m))
    """

    # Initialise an empty Grammaian, the parameter space and calculate the size of the
    # parameter space (total number of points)
    Grammian = np.array([])
    paramSpace = []
    totalNumber = 1
    for i in params:
        paramSpace += [ np.linspace(i[0], i[1], i[2]) ]
        totalNumber = totalNumber * i[2]
    
    # Seed choice (arbitrary)
    RB = np.zeros(len(paramSpace))
    for i in range(RB.shape[0]):
        # FIXME randomize this bit.
        RB[i] = paramSpace[i][0]

    # Initiate the greedy algorithm
    sigma = np.append(1)

    # The parameter space is large, so the computation will be expensive
    while sigma[-1] > epsilon:
        i = RB.shape[0]
        sigma = np.append(0)
        Grammian = constructGrammian (RB, matrix, Grammian)

        # Construct the Gram matrix inverse
        G_inv = np.inverse(Grammian)

        # FIXME I need to check the rest of this function thoroughly
        # Stupidly traverse the entire parameter space
        # NOTE: We could have MCMC or a simple MC method as well
        # Can we edit the ranges where we are searching by discarding regions in
        # parameter space when we find the vectors? (Suggestion by Priscilla)
        for j in range(totalNumber):
            # the params_trial does not have to be zeroed away
            # This is a temporary variable
            translateIt = j

            # Construct the test vector of the parameters
            for k in range(len(paramSpace)):
                n = translateIt % params[k][2]
                params_trial[k] = paramSpace[k][n]
                translateIt = translateIt // params[k][2]

            # Calculate the projection
            projection = 0
            for k in range(RB.shape[0]):
                for l in range(RB.shape[0])):
                    projection += G_inv[k,l] * innerProduct(
                                sm.dataGeneration(schedule, params_trial, pulsars)
                                sm.dataGeneration(schedule, RB[l], pulsars) 
                                matrix) * sm.dataGeneration(schedule, RB[k], pulsars)

            norm = sm.dataGeneration(schedule, params_trial, pulsars) - projection
            norm *= np.conjugate(norm)

            if sigmaTrial > sigma[i]:
                sigma[i] = sigmaTrial
                params_i = params_trial

        # Add the lambda_i, which was found by maximizing
        RB = np.append(RB, params_i)

    return RB

def covarianceMatrix (pulsars, schedule, GWB = false, WhiteNoise = true, RedNoise = false, PowerLaw = false):
    # Construct a zero matrix
    matrix = np.zeros (vec1.shape, vec1.shape)
    cdef int a, b, N
    cdef double ta, tb, gamma, gamma_a, tau, f_L, freqerror

    # Set gamma to 7/3 as we assume the the GWB comes mainly from SMBHBs
    gamma = 7/3
    # The lower frequency cut off. It should be much larger than the observation time
    freqerror = 1e-5
    f_L = freqerror / pulsars[1,-1]

    # Truncate the series in the power law and the GWB noise after N iterations
    N = 1e4

    # Let the GWB amplitude be small
    A_GWB = 0.01

    if GWB:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                a, ta = schedule[:,i]
                b, tb = schedule[:,j]
                tau = 2*np.pi * (ta - tb)
                # From the pulsars we need only the unit vectors and the noise
                # magnitude, we could probably do it here, as it would make the code
                # look better
                C = 1/2 * (1 - np.dot( pulsars.getUnitVector(a), pulsars.getUnitVector(b)))
                matrix[i,j] += sm.covarianceMatrixMemberGWB (i, j, a, b, A, f_L, gamma, tau,
                        N, C)

    if WhiteNoise:
        for i in range(matrix.shape[0]):
            a = schedule[0,i]
            A_WN = pulsars.getWhiteNoise(a)
            matrix[i,i] += sm.covarianceMatrixMemberWN (i, i, a, a, A_WN)

    # FIXME Implement the methods
    # if RedNoise
    # if PowerLaw

    return matrix


def innerProduct (vec1, vec2, matrix):
    """
    This is the vector product, which is the main bottleneck in the calculation.

    The function needs to be given two vectors and a matrix and the dimensionalities
    should match.
    """
    # The return value
    r = 0

    if matrix.shape[0] == vec2.shape[0] and matrix.shape[1] == vec1.shape[0]:
        # Do a vector multiplication (this method should not use too much memory)
        for i in ranger(vec2.shape[0]):
            r += vec1[i] * np.dot(matrix[:,i],vec2)
    else:
        print("Error!!!!")

    return r

def constructGrammian (set, matrix, G=np.array([])):
    """
    Here we construct the Grammian matrix and we do it in a clever way. We take the
    previous matrix and extend it, because when we increase the number of basis in our
    set, we need to change only a very small part of the matrix.

    Also we are using the fact that the Grammian matrix is symmetric.

    Also, if the given Grammian is empty, but RB is larger than one member, we will
    calculate the Grammian from scratch.
    """
    tmp = G
    G = np.zeros(set.shape[0], set.shape[0])
    calculateAll = (G.shape[0] != 1 and tmp.shape[0] == 0)
    if not calculateAll and tmp.shape[0] == 0:
        G[:-1,:-1] = tmp

    if not calculateAll:
        # Generate one of the templates
        vec1 = sm.dataGeneration (schedule, set[-1], pulsars)
        G[-1,-1] = innerProduct(template1, template1)

        # Use the fact that Grammian is symmetric
        for i in range(G.shape[0] - 1):
            vec2 = sm.dataGeneration(schedule, set[i], pulsars)
            G[i,-1] = innerProduct(vec1, vec2, matrix)
            G[-1,i] = G[i,-1]
    else:
        # Use the fact that Grammian is symmetric
        for i in range(G.shape[0]):
            # Generate one of the templates
            vec1 = sm.dataGeneration (schedule, set[i], pulsars)
            for j in range(G.shape[0]):
                vec2 = sm.dataGeneration (schedule, set[j], pulsars)
                G[i,j] = innerProduct(vec1, vec2, matrix)
                if i != j:
                    G[j,i] = G[i,j]

    return G
