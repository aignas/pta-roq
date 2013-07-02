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
from libc.math cimport sin, cos, sqrt, exp, log

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

        # Three kinds of noises. And the actual realisations of them are not correlated
        self.__Noise = np.random.rand(N,3)

        # Transform the random number distributions to span the ranges defined by the
        # user. Also scale the noises so that we can define different ratios of white,
        # red and power law noise.
        for i in range(3):
            self.__R[:,i] = self.__R[:,i] * abs(ranges[i,0] - ranges[i,1]) + ranges[i].min()
            self.__Noise[:,i] *= noise[i]

    def getNoise (self, int idx):
        return self.__Noise[idx,:]

    def getWhiteNoise (self, int idx):
        return self.__Noise[idx,0]

    def getRedNoise (self, int idx):
        return self.__Noise[idx,1]

    def getPowLawNoise (self, int idx):
        return self.__Noise[idx,2]

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
def pulsarTerm (double t, double M, double D, double iota, double Phi0, double psi,
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
def noise (double t, np.ndarray variance):
    cdef double white, red, powerLaw

    # Implement the white  noise
    white = np.random.rand() * sqrt(variance[0])

    # FIXME Implement the red and power law noises as well as shown in the van Haasteren
    # et al. paper from 2009
    # I found this: http://stackoverflow.com/questions/918736 , but do not know yet how
    # to apply it.

    # Turn these two types of noise of
    red = 0
    powerLaw = 0

    # Add all the noise contributions
    return white + red + powerLaw

# Define the residual as a function of parameters
def individualSource (double t, np.ndarray params, double L, np.ndarray u_p):
    cdef double M, D, iota, Phi0, psi, theta, phi, omega0, zeta, p
    cdef np.ndarray a, A

    # Unpack all the parameters:
    M, D, iota, Phi0, psi, theta, phi, omega0 = params.tolist()

    zeta = M**(5/3)/D
    a = amplitude(zeta, iota, phi, psi)
    A = basis (omega0, M, theta, phi, t, u_p)
    p = pulsarTerm (t, M, D, iota, Phi0, psi, theta, phi, omega0, L, u_p)
    return np.dot(a,A) + p

# Define the residual as a function of parameters
def residual (double t, np.ndarray sources, double L, np.ndarray u_p, np.ndarray variance):
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
def dataGeneration (np.ndarray schedule, np.ndarray sources, pulsars, 
        addNoise=True, addGWB=False):
    cdef double t, L
    cdef np.ndarray u_p, a = np.array([]), noise = np.zeros(3)
    cdef int N = pulsars.getNumber()

    # Copy the structure of the array for the time log

    # Start collecting the data
    for i in range(schedule.shape[1]):

        # Set some temporary values
        L = pulsars.getLength(schedule[0,i])
        u_p = pulsars.getUnitVector(schedule[0,i])

        if addNoise:
            noise = pulsars.getNoise(schedule[0,i])

        a = np.append(a, residual(schedule[1,i], sources, L, u_p, noise))

    return a

# Calculate the GWB terms in the covariance matrix members.
def covarianceMatrixMemberGWB (i, j, a, b, A, f, gamma, tau, N, C):
    cdef double alpha, sum, sum_member

    alpha = 3/2 * C * log(C) - C/4 + 1/2
    # a simple delta function implementation
    if a == b:
        alpha += 1/2

    # Here I use slightly more memory by storring each_member before summing it, but
    # this way I do not have to calculate horrible factorials and it should speed things
    # up a bit
    # This function calculates N terms and then truncates the series
    sum_member = - 1 / (1 + gamma)
    sum = sum_member
    for i in range(N):
        sum_member *= - (f * tau)**2 / ((2*i + 1) * (2*i + 2)) * (2*i - 1 - gamma) \
                      / (2*i + 1 - gamma)
        sum += sum_member

    return A**2 * alpha / ((2 * np.pi)**2 * f**(1 + gamma)) \
           / (np.gamma(-1 - gamma) * sin(-np.pi * gamma / 2) * (f*tau)**(1 + gamma) - sum)

# Calculate white noise terms
def covarianceMatrixMemberWN (i,j,a,b,N):
    # r is the return value and N should be the noise amplitude
    cdef double r = 0

    if a==b and i==j:
        r = N**2

    return r

# Calculate red noise Lorentzian terms
def covarianceMatrixMemberLor (i,j,a,b,N,f,tau):
    # r is the return value and N should be the noise amplitude
    cdef double r = 0
    
    if a==b:
        r = N**2 * exp(-f*tau)

    return r

# Calculate the power law spectral noise
def covarianceMatrixMemberPowLaw (i, j, a, b, A, f, tau, gamma, N):
    cdef double r = 0

    if a==b:
        # Here I use a similar technique to the one explained in
        # covarianceMatrixMemberGWB
        sum_member = 1 / (1 - gamma)
        sum = sum_member
        for i in range(N):
            sum_member *= - (f * tau)**2 / ((2*i + 1) * (2*i + 2)) * (2*i + 1 - gamma) \
                          / (2*i + 3 - gamma)
            sum += sum_member

        r =  A**2 / (f**(gamma - 1)) \
             * (np.gamma(1 - gamma) * sin(np.pi * gamma /2) * (f * tau)**( gamma -1) - sum)
    
    # Return the computed value
    return r
