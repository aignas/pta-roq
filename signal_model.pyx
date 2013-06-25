#! /usr/bin/python2.7

# Have some valuable things from Python 3
from __future__ import print_function, division

# The following was copied from a Cython manual
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.int
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_t
# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.
from libc.math cimport sin, cos, sqrt

# The following will code up the equations (16) to (19) from the Ellis et al. (see
# README)

# This is a class for unit vectors. I have a feeling, that the p vector might be not
# necessary here
class UnitVectors:
    def __init__ (self, double theta, double phi):
        # Shorthands for sines and cosines
        self.theta = theta
        self.phi = phi

    def Omega (self):
        # define a unit vector \hat{\Omega} as shown in the (4)
        return np.array([
        -sin(self.theta)*cos(self.theta),
            -sin(self.theta)*sin(self.phi),
            -cos(self.theta)])

    def m (self):
        # define a unit vector \hat{m} as shown in the (5)
        return np.array([
            -sin(self.phi),
            cos(self.phi),
            0])

    def n (self):
        # define a unit vector \hat{n} as shown in the (6)
        return np.array([
            -cos(self.theta)*cos(self.theta),
            -cos(self.theta)*sin(self.phi),
            sin(self.theta)])

# This will be a Pulsar Array data structure and we access the data in it the lazy
# way.
class PulsarGrid:
    def __init__ (self, int N, ranges, noise):
        # These wil be polar coordinates of the pulsars
        self.__R = np.random.rand(N,3)
        self.__Noise = np.random.rand(N)*noise

        # Transform the random number distributions to span the ranges defined by the
        # user
        # FIXME Is there a nicer way to write it?
        for i in range(3):
            self.__R[:,i] = self.__R[:,i] * abs(ranges[i,0] - ranges[i,1]) + ranges[i].min()

    def getNoise (self, int idx):
        return self.__Noise[idx]

    def getAngles (self):
        return self.__R[:,1:3]

    def getUnitVector (self, int idx):
        p = np.zeros(3)
        p[0] =  sin (self.__R[idx,1])
        p[1] =  p[0]
        p[2] =  cos (self.__R[idx,1])
        p[0] *= cos (self.__R[idx,2])
        p[1] *= sin (self.__R[idx,2])
        return p

    def getLength (self, int idx):
        return self.__R[idx,0]

    def getNumber (self):
        return self.__R.shape[0]

# define antenna pattern functions as in (9) and (10)
# The F[0] term is Fplus and F[1] is Fcross
def antennaPattern (double theta, double phi, np.ndarray u_p):
    # Initialise unit vector according to the given parameters
    v = UnitVectors(theta, phi)

    cdef double d1, d2, d3
    cdef np.ndarray F

    d1 = np.dot(v.m(),u_p)
    d2 = np.dot(v.n(),u_p)
    d3 = np.dot(v.Omega(),u_p)

    F = np.zeros(2)
    F[0] = 0.5 * ( d1 + d2 ) * ( d1 - d2 ) / ( 1 + d3 )
    F[1] = ( d1 * d2 ) / ( 1 + d3 )

    return F

# Define some required functions.
# FIXME we can probably do some approximations here and there (See the Ellis et al.
# paper p3).
def omega(double t, double omega0, double M):
    return (omega0**(-8/3) - 256/5 * M**(5/3) * t)**(-3/8)

def Phi (double t, double omega0, double M):
    # FIXME the Cornell et al. paper has a constant term here. I believe, that the
    #       constant term here should be as well, but I am wondering if the constant
    #       term is the same as \Phi_0
    return 1/(32 * M**(5/3)) * ( omega0**(-5/3) - omega(t, omega0, M)**(-5/3))

# define the GW contributions to the timing residuals as in (12) and (13)
# The first term is the plus, the second term is the cross as in the F function
def gWContribution (double omega0, double M, double t, double iota, double psi,
        double zeta, double Phi0):

    cdef double a, b
    cdef np.ndarray gw

    a = sin( 2 * (Phi(t, omega0, M) - Phi0) ) * (1 + cos(iota)**2)
    b = cos( 2 * (Phi(t, omega0, M) - Phi0) ) * cos(iota)
    gw = np.zeros(2)

    gw[0] = zeta / omega(t, omega0, M)**(1/3) * ( - a * cos(2*psi) - 2 * sin(2*psi) )
    gw[1] = zeta / omega(t, omega0, M)**(1/3) * ( - a * sin(2*psi) - 2 * cos(2*psi) )

    return gw


# define the coefficients amplitudes as shown in the (18)
def amplitude(double zeta, double iota, double phi, double psi):
    cdef np.ndarray a

    a = np.zeros(4)
    a[0] = zeta * ( (1 + cos(iota)**2) * cos (phi) * cos (2*psi) \
                + 2 * cos(iota) * sin (phi) * sin (2*psi) )

    a[1] = - zeta * ( (1 + cos(iota)**2) * sin (phi) * cos (2*psi) \
                + 2 * cos(iota) * cos (phi) * sin (2*psi) )

    a[2] = zeta * ( (1 + cos(iota)**2) * cos (phi) * sin (2*psi) \
                + 2 * cos(iota) * sin (phi) * cos (2*psi) )

    a[3] = - zeta * ( (1 + cos(iota)**2 ) * sin (phi) * sin (2*psi) \
                + 2 * cos(iota) * cos (phi) * cos (2*psi) )

    return a


# Define the time dependent basis functions as shown in the equation (19)
def basis (double omega0, double M, double theta, double phi, double t, np.ndarray u_p):
    cdef np.ndarray A, F

    A = np.zeros(4)
    F = antennaPattern(theta, phi, u_p)
    A[0] = F[0] * omega(t, omega0, M)**(-1/3) * sin (2 * Phi(t, omega0, M))
    A[1] = F[0] * omega(t, omega0, M)**(-1/3) * cos (2 * Phi(t, omega0, M))
    A[2] = F[1] * omega(t, omega0, M)**(-1/3) * sin (2 * Phi(t, omega0, M))
    A[3] = F[1] * omega(t, omega0, M)**(-1/3) * cos (2 * Phi(t, omega0, M))

    return A

# Define the pulsar term as in the eq (17)
def pulsar (double t, double M, double D, double iota, double Phi0, double psi,
        double theta, double phi, double omega0, double L, np.ndarray u_p):
    cdef np.ndarray F, s
    cdef double tp, zeta

    v = UnitVectors(theta, phi)
    tp = t - L * (1 + np.dot(v.Omega(), u_p))
    zeta = M**(5/3)/D

    F = antennaPattern(theta, phi, u_p)
    s = gWContribution (omega0, M, tp, iota, psi, zeta, Phi0)
    return np.dot(F,s)

# Define the noise term
def noise (double t, double var):
    # Add some gaussian noise which deppends on each pulsar
    return np.random.randn() * sqrt(var)

# Define the residual as a function of parameters
def individualSource (double t, np.ndarray params, double L, np.ndarray u_p):
    cdef double M, D, iota, Phi0, psi, theta, phi, omega0, zeta, p
    cdef np.ndarray a, A

    # Unpack all the parameters:
    M, D, iota, Phi0, psi, theta, phi, omega0 = params.tolist()

    zeta = M**(5/3)/D
    a = amplitude(zeta, iota, phi, psi)
    A = basis (omega0, M, theta, phi, t, u_p)
    p = pulsar (t, M, D, iota, Phi0, psi, theta, phi, omega0, L, u_p)
    return np.dot(a,A) + p

# Define the residual as a function of parameters
def residual (double t, np.ndarray sources, double L, np.ndarray u_p, double variance):
    cdef double Signal = 0

    # Add contributions to the signal from all GW sources
    for i in range(sources.shape[0]):
        Signal += individualSource(t, sources[i], L, u_p)

    # Add noise to the signal
    Signal += noise (t, variance)

    return Signal

def genSchedule (np.ndarray schedule, double t_final, double dt_min, double dt_max):
    cdef double t
    cdef np.ndarray u_p
    cdef np.ndarray dates_out = np.array([]), index_out = np.array([])
    cdef int N = schedule.shape[0]
    collectData = True

    # Copy the structure of the array for the time log
    dates = [np.array([])]*N

    # Start collecting the data
    while collectData:
        collectData = False

        for i in range(N):
            t = schedule[i]

            if t > t_final:
                continue
            else:
                # The data is being collected, add to the log
                dates[i] = np.append(dates[i], t)
                # Do not stop generating data if there is at least one pulsar, which had
                # not "gone" in to the future
                collectData = True

            # Update a schedule for this pulsar. We randomise the next measurement a bit
            schedule[i] += np.random.rand() * abs(dt_max - dt_min) + dt_min

    # FIXME Spit out the data in a format we need
    # Contract the data into a one vector
    for i in range(N):
        dates_out = np.append(dates_out, dates[i])
        index_out = np.append(index_out, [i]*dates[i].shape[0])

    return np.append([index_out], [dates_out], axis=0)

# Generate the actual data.
def dataGeneration (np.ndarray schedule, np.ndarray sources, pulsars, addNoise = True):
    cdef double t, L, noise
    cdef np.ndarray u_p
    cdef np.ndarray a = np.array([])
    cdef int N = pulsars.getNumber()

    # Copy the structure of the array for the time log

    # Start collecting the data
    for i in range(schedule.shape[1]):

        # Set some temporary values
        L = pulsars.getLength(schedule[0,i])
        u_p = pulsars.getUnitVector(schedule[0,i])
        if addNoise:
            noise = pulsars.getNoise(schedule[0,i])
        else:
            noise = 0
        a = np.append(a, residual(schedule[1,i], sources, L, u_p, noise))

    return a
