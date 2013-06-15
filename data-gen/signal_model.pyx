#! /usr/bin/python2.7

import numpy as np
from libc.math cimport sin, cos

# The following will code up the equations (16) to (19) from the Ellis et al. (see
# README)

# This is a class for unit vectors. I have a feeling, that the p vector might be not
# necessary here
class UnitVectors:
    def __init__ (self, theta, phi):
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
    def __init__ (self, N):
        # These wil be polar coordinates of the pulsars
        self.__R = np.random.rand(N*3)
        self.__R = self.__R.reshape(N,3)

        # FIXME the 10 here needs to be a general length scale 5 as well
        self.__R[0,:] = self.__R[0,:] * 10 + 5
        # Angular Parts range from -pi,pi or 0,pi
        self.__R[1,:] = (self.__R[1,:] - 0.5) * 2*np.pi
        self.__R[2,:] = self.__R[2,:] * np.pi

    def getAngles (self, idx):
        return self.__angles[:,idx]

    def getUnitVector (self, idx):
        p = np.array([sin (self.__angles[0,idx]), 0, 0 ])
        p[1] = p[0]
        p[2] = np.cos (self.__angles[0,idx])
        p[0] *= np.cos (self.__angles[1,idx])
        p[1] *= np.sin (self.__angles[1,idx])
        return p

    def getLength (self, idx):
        return self.__Length[idx]

# define antenna pattern functions as in (9) and (10)
# The F[0] term is Fplus and F[1] is Fcross
def antennaPattern (theta, phi):
    # Initialise unit vector according to the given parameters
    v = UnitVectors(theta, phi)

    F = np.array([0,0])
    F[0] = 0.5 * ( np.dot(v.m,v.p)**2 - np.dot(v.n,v.p)**2 )/ \
              ( 1 + np.dot(v.Omega,v.p))
    F[1] = ( np.dot(v.m,v.p) * np.dot(v.n,v.p) )/( 1 + np.dot(v.Omega,v.p))
    
    return F

# Define some required functions.
# FIXME we can probably do some approximations here and there (See the Ellis et al.
# paper p3).
def omega(t, omega0, M):
    return (omega0**(-8./3) - 256./5 * M**(5./3) * t)**(-3./8)

def Phi (t, omega0, M):
    # FIXME the Cornell et al. paper has a constant term here. I believe, that the
    #       constant term here should be as well, but I am wondering if the constant
    #       term is the same as \Phi_0
    return 1./(32 * M**(5./3)) * ( omega0**(-5./3) - omega(t)**(-5./3))

# define the GW contributions to the timing residuals as in (12) and (13)
# The first term is the plus, the second term is the cross as in the F function
def gWContribution (t, iota, psi, zeta, Phi0):
    a = sin( 2 * (Phi(t) - Phi0) ) * (1 + cos(iota)**2)
    b = cos( 2 * (Phi(t) - Phi0) ) * cos(iota)
    gw = np.array([0,0])

    gw[0] = zeta / omega(t)**(1./3) * ( - a * cos(2*psi) - 2 * sin(2*psi) )
    gw[1] = zeta / omega(t)**(1./3) * ( - a * sin(2*psi) - 2 * cos(2*psi) )

    return gw


# define the coefficients amplitudes as shown in the (18)
def amplitude(zeta, iota, phi, psi):
    a = np.array([ 0,0,0,0 ])
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
def basis (theta, phi, t):
    A = np.arraf([0,0,0,0])
    F = antennaPattern(theta, phi, t)
    A[0] = F[0] * omega(t)**(-1./3) * sin (2 * Phi(t))
    A[1] = F[0] * omega(t)**(-1./3) * cos (2 * Phi(t))
    A[2] = F[1] * omega(t)**(-1./3) * sin (2 * Phi(t))
    A[3] = F[1] * omega(t)**(-1./3) * cos (2 * Phi(t))

    return A

# Define the pulsar term as in the eq (17)
def pulsar (t, zeta, iota, Phi0, psi, theta, phi, omega0, L):
    v = UnitVectors(theta, phi)
    tp = t - L * (1 + np.dot(v.Omega, v.p))

    F = antennaPattern(theta, phi, t)
    s = gWContribution (tp, iota, psi, zeta, Phi0)
    return np.dot(F,s)

# Define the noise term
def noise (t):
    #define some Gaussian noise in the amplitude.
    var = 0.00003
    return np.exp ( - np.random.rand()**2 / (2*var))

# Define the residual as a function of parameters
def residual (t, zeta, iota, Phi0, psi, theta, phi, omega0, L):
    a = amplitude(zeta, iota, phi, psi)
    A = basis (theta, phi, t)
    p = pulsar (t, zeta, iota, Phi0, psi, theta, phi, omega0, L)
    n = noise (t)
    return np.dot(a,A) + p + n
