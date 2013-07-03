#! /usr/bin/python2.7

from __future__ import print_function, division

# Check out signal_model.pyx for explanation of the following lines
import numpy as np
cimport numpy as np

from libc.math cimport sin, cos, sqrt

import pyximport; pyximport.install()
import signal_model as sm

from cython.parallel import parallel, prange

def ROMGreedy (np.ndarray schedule, pulsars, np.ndarray[double, ndim=2] params, \
        np.ndarray[double, ndim=2] matrix, double epsilon):
    """This will generate a reduced basis in the given parameter space.

    Args:

        schedule: An array which contains all the ids of the pulsars and when they
        were/going to be measured.

        pulsars: The pulsar array, which contains all the properties of the pulsars

        params: A parameter space value. At the moment it is assumed to be of the
            following structure:
                ((param^1_min, param^1_max, param^1_step),
                 (param^2_min, param^2_max, param^2_step),
                 ...
                 (param^n_min, param^n_max, param^n_step))

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
    cdef np.ndarray[double, ndim=2] Grammian, Grammian_inv, RB, RB_tmp, \
            projectionCoeffs, projectionCoeffs_tmp
    cdef np.ndarray[double] sigma, params_trial, sigma_trial, paramSpace, params_tmp, \
            projection_tmp, projection
    cdef np.ndarray paramDimensions
    cdef int i

    # Initialize empty arrays (or allocate memory):
    Grammian = np.array([[]])
    paramDimensions = np.array([], dtype='i4')
    paramSpace = np.array([])

    for i in xrange(params.shape[0]):
        paramSpace = np.append(paramSpace, np.linspace(params[i,0], params[i,1], params[i,2]))
        paramDimensions = np.append(paramDimensions, int(params[i,2]))

    totalNumber = paramDimensions.prod()

    sigma_trial = np.zeros(totalNumber)
    projectionCoeffs = np.zeros((totalNumber,1))
    
    # Seed choice (arbitrary). We just randomize the first choice. We also select the
    # first error arbitrary. This is just to have the same dimensions of two arrays
    sigma = np.array([1.], dtype=np.double)

    # Randomize the first basis
    RB = np.zeros((params.shape[0],1))
    params_tmp = ID2param(np.random.randint(totalNumber), paramSpace, paramDimensions)
    for i in xrange(RB.shape[0]):
        RB[i,0] = params_tmp[i]

    # The parameter space is large, so the computation will be expensive
    while sigma[-1] > epsilon:
        # Construct the Gram matrix and its inverse
        Grammian, Grammian_inv = constructGrammian (schedule, pulsars, RB, matrix, Grammian)

        # Stupidly traverse the entire parameter space
        # NOTE: We could have MCMC or a simple MC method as well
        # Can we edit the ranges where we are searching by discarding regions in
        # parameter space when we find the vectors? (Suggestion by Priscilla)
        for i in xrange(totalNumber):
            # the params_trial does not have to be zeroed away
            # This is a temporary variable
            params_trial = ID2param(i, paramSpace, paramDimensions)

            # Calculate the projection
            projection_tmp = projectOnRBcoeff(schedule, pulsars, RB, matrix,
                    params_trial, Grammian_inv, projectionCoeffs[i])
            projectionCoeffs[i,-1] = projection_tmp[-1]

            projection = projectOnRB(schedule, pulsars, RB, projection_tmp)

            # Calculate modulus in a clever way (Or maybe not so clever, but I believe,
            # that it will use slightly less memory comparing to a oneliner)
            difference = sm.dataGeneration(schedule, params_trial, pulsars) - projection
            sigma_trial[i] = np.dot(difference,difference)

        sigma = np.append(sigma, sigma_trial.max())
        N_i = sigma_trial.argmax()

        # Add the lambda_i, which was found by maximizing the error
        params_tmp = ID2param(N_i, paramSpace, paramDimensions)
        RB_tmp = RB
        RB = np.zeros((RB_tmp.shape[0], RB_tmp.shape[1] + 1 ))
        for i in xrange(RB_tmp.shape[0]):
            RB[i,RB_tmp.shape[1]] = params_tmp[i]
            for j in xrange(RB_tmp.shape[1]):
                RB[i,j] = RB_tmp[i,j]

        # extend the projectionCoeffs matrix
        projectionCoeffs_tmp = projectionCoeffs
        projectionCoeffs = np.zeros((projectionCoeffs.shape[0], projectionCoeffs.shape[1] + 1))
        for i in xrange(projectionCoeffs_tmp.shape[0]):
            for j in xrange(projectionCoeffs_tmp.shape[1]):
                projectionCoeffs[i,j] = projectionCoeffs_tmp[i,j]

        print("Number of RB and the error: " + str(sigma.shape[0]) + " " + str(sigma[-1]))

    return RB

def projectOnRBcoeff (np.ndarray schedule, pulsars, np.ndarray[double, ndim=2] RB, 
        np.ndarray[double, ndim=2] matrix, np.ndarray[double] params_trial,
        np.ndarray[double, ndim=2] G_inv,  np.ndarray[double] coefs):
    # Perform various checks whether the given parameters are valid.
    if G_inv.shape[0] != RB.shape[1]:
        print("ERROR: The given grammian doesn't have the required dimensions")

    # Initialize some variables
    cdef int i = coefs.shape[0] - 1, j

    for j in xrange(RB.shape[1]):
        coefs[i] += G_inv[i,j] * innerProduct(
                sm.dataGeneration(schedule, params_trial, pulsars),
                sm.dataGeneration(schedule, RB[:,j], pulsars),
                matrix)

    return coefs

# FIXME this function doesn't work
def projectOnRB (np.ndarray schedule, pulsars, np.ndarray[double, ndim=2] RB,
        np.ndarray[double, ndim=1] coefs):
    # Initialize some variables
    cdef np.ndarray[double] p = coefs[0] * sm.dataGeneration(schedule, RB[:,0], pulsars)
    cdef int i, j

    for i in xrange(1,RB.shape[1]):
        p += coefs[i] * sm.dataGeneration(schedule, RB[:,i], pulsars)

    return p

def covarianceMatrix (pulsars, np.ndarray schedule, GWB=False, 
        WhiteNoise=True, RedNoise=False, PowerLaw=False):
    # Construct a zero matrix
    cdef np.ndarray[double, ndim=2] matrix = np.zeros((schedule.shape[1], schedule.shape[1])), matrix_inv
    cdef int a, b, N
    cdef double ta, tb, gamma, gamma_a, tau, f_L, freqerror
    cdef int i, j

    # Set gamma to 7/3 as we assume the the GWB comes mainly from SMBHBs
    gamma = 7/3
    # The lower frequency cut off. It should be much larger than the observation time
    freqerror = 1e-5
    f_L = freqerror / schedule[1,-1]

    # Truncate the series in the power law and the GWB noise after N iterations
    N = 10000

    if GWB:
        # Let the GWB amplitude be small
        A_GWB = 0.01
        for i in xrange(matrix.shape[0]):
            for j in xrange(matrix.shape[1]):
                a, ta = schedule[:,i]
                b, tb = schedule[:,j]
                tau = 2*np.pi * (ta - tb)
                # From the pulsars we need only the unit vectors and the noise
                # magnitude, we could probably do it here, as it would make the code
                # look better
                C = 1/2 * (1 - np.dot( pulsars.getUnitVector(a), pulsars.getUnitVector(b)))
                matrix[i,j] += sm.covarianceMatrixMemberGWB (i, j, a, b, A_GWB, f_L, gamma, tau,
                        N, C)

    if WhiteNoise:
        for i in xrange(matrix.shape[0]):
            a = schedule[0,i]
            A_WN = pulsars.getWhiteNoise(a)
            matrix[i,i] += sm.covarianceMatrixMemberWN (i, i, a, a, A_WN)

    # FIXME Implement the methods
    # if RedNoise
    # if PowerLaw

    # Invert the matrix
    matrix_inv = np.linalg.inv(matrix)

    return matrix, matrix_inv


def innerProduct (np.ndarray vec1, np.ndarray vec2, np.ndarray matrix):
    """
    This is the vector product, which is the main bottleneck in the calculation.

    The function needs to be given two vectors and a matrix and the dimensionalities
    should match.
    """
    # Initialise the return value
    cdef double r

    if matrix.shape[0] == vec2.shape[0] and matrix.shape[1] == vec1.shape[0]:
        # Do a vector multiplication (this method should not use too much memory)
        r = np.dot(vec1, np.dot(matrix, vec2))

    else:
        print("Error!!! The inner product is not working because of the missmatch of dimensions")

    return r

def constructGrammian (np.ndarray schedule, pulsars, np.ndarray[double, ndim=2] set,
        np.ndarray[double, ndim=2] matrix, np.ndarray[double, ndim=2] tmp=np.array([[]])):
    """
    Here we construct the Grammian matrix and we do it in a clever way. We take the
    previous matrix and extend it, because when we increase the number of basis in our
    set, we need to change only a very small part of the matrix.

    Also we are using the fact that the Grammian matrix is symmetric.

    Also, if the given Grammian is empty, but RB is larger than one member, we will
    calculate the Grammian from scratch.

    Args:
        set: The basis set to use for construction

        matrix: The matrix to use for the inner product

        G: The Grammian matrix which we want to extend, because the set was only
        extended from the previous calculation

    Returns:

        New G and the inverse of it.
    """
    cdef np.ndarray[double, ndim=2] G, G_inv
    cdef np.ndarray[double] vec1, vec2, params1, params2
    cdef int i, j, k, size = tmp.shape[1], paramSize = set.shape[0]

    G = np.zeros((set.shape[1], set.shape[1]))

    # Check if the given G dimensionality is one lower than the set and if the given set
    # consists more than one vector

    if (size == 0):
        print("Constructing the Gram matrix")
        # Use the fact that Grammian is symmetric
        for i in xrange(G.shape[0]):

            # Construct a parameter vector
            params1 = np.zeros(paramSize)
            for j in xrange(paramSize):
                params1[j] = set[j,i]

            # Generate one of the templates
            vec1 = sm.dataGeneration (schedule, params1, pulsars)

            G[i,i] = innerProduct(vec1, vec1, matrix)
            for j in xrange(i+1, G.shape[0]):

                # Construct a parameter vector
                params2 = np.zeros(paramSize)
                for k in xrange(paramSize):
                    params2[k] = set[k,j]

                vec2 = sm.dataGeneration (schedule, params2, pulsars)
                G[i,j] = innerProduct(vec1, vec2, matrix)
                G[j,i] = G[i,j]
    else:
        print("Extending the Gram matrix")
        for i in xrange(size):
            G[i,i] = tmp[i,i]
            for j in xrange(i + 1, size - 1):
                G[i,j] = tmp[i,j]
                G[j,i] = tmp[j,i]

        # Construct a parameter vector
        params1 = np.zeros(paramSize)
        for j in xrange(paramSize):
            params1[j] = set[j,size]

        # Generate one of the templates
        vec1 = sm.dataGeneration (schedule, params1, pulsars)
        G[size,size] = innerProduct(vec1, vec1, matrix)

        # Use the fact that Grammian is symmetric
        for i in xrange(size):

            # Construct a parameter vector
            params2 = np.zeros(paramSize)
            for j in xrange(paramSize):
                params2[j] = set[j,size]

            vec2 = sm.dataGeneration(schedule, params2, pulsars)
            G[i][size] = innerProduct(vec1, vec2, matrix)
            G[size][i] = G[i][size]

    print("Done")
    # Calculate the inverse of a matrix
    G_inv = np.linalg.inv(G)

    return G, G_inv

def ID2param (int idx, np.ndarray[double] space, np.ndarray[long] dim):
    # Initialize a return container
    cdef np.ndarray r = np.zeros(dim.shape[0])
    cdef int n, i

    # Construct the test vector of the parameters
    for i in xrange(r.shape[0]):
        n = dim[0:i].sum() + idx % dim[i]
        r[i] = space[n]
        idx = idx // dim[i]

    return r
