#! /usr/bin/python2.7

import numpy as np
from libc.math cimport sin, cos

# The following will code up the equations (16) to (19) from the Ellis et al. (see
# README)

# This is a class for unit vectors. I have a feeling, that the p vector might be not
# necessary here
class unit_vectors:
    def __init__ (self, theta, phi, theta_p, phi_p):
        # Shorthands for sines and cosines
        st = sin(theta)
        sp = sin(phi)
        ct = cos(theta)
        cp = cos(phi)

        # define a unit vector \hat{\Omega} as shown in the (4)
        self.Omega = np.array([-st*cp, -st*sp, -ct])
        # define a unit vector \hat{m} as shown in the (5)
        self.m = np.array([-sp, cp, 0])
        # define a unit vector \hat{n} as shown in the (6)
        self.n = np.array([-ct*cp, -ct*sp, st])

        # unit vector pointing from the earth to the pulsar
        # FIXME check whether this term is correct
        self.p = np.array([
            sin(theta_p)*cos(phi_p), 
            sin(theta_p)*sin(phi_p),
            cos(theta_p)
            ])

# define antenna pattern functions as in (9) and (10)
# The F[0] term is Fplus and F[1] is Fcross
def AntennaPattern (theta, phi, theta_p, phi_p):
    # Initialise unit vector according to the given parameters
    v = unit_vectors(theta, phi, theta_p, phi_p)

    F = []
    F += [ 0.5 * ( np.dot(v.m,v.p)**2 - np.dot(v.n,v.p)**2 )/ \
            ( 1 + np.dot(v.Omega,v.p)) ]
    F += [ ( np.dot(v.m,v.p) * np.dot(v.n,v.p) )/( 1 + np.dot(v.Omega,v.p)) ]
    
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
def s (t, iota, psi, zeta, Phi0):
    a = sin( 2 * (Phi(t) - Phi0) ) * (1 + cos(iota)**2)
    b = cos( 2 * (Phi(t) - Phi0) ) * cos(iota)
    gw = []

    gw += [ zeta / omega(t)**(1./3) * ( - a * cos(2*psi) - 2 * sin(2*psi) ) ]
    gw += [ zeta / omega(t)**(1./3) * ( - a * sin(2*psi) - 2 * cos(2*psi) ) ]

    return gw


# define the coefficients amplitudes as shown in the (18)
def amplitude(zeta, iota, phi, psi):
    a = []
    a += [ zeta * ( (1 + cos(iota)**2) * cos (phi) * cos (2*psi) \
                + 2 * cos(iota) * sin (phi) * sin (2*psi) ) ]

    a += [ - zeta * ( (1 + cos(iota)**2) * sin (phi) * cos (2*psi) \
                + 2 * cos(iota) * cos (phi) * sin (2*psi) ) ]

    a += [ zeta * ( (1 + cos(iota)**2) * cos (phi) * sin (2*psi) \
                + 2 * cos(iota) * sin (phi) * cos (2*psi) ) ]

    a += [ - zeta * ( (1 + cos(iota)**2 ) * sin (phi) * sin (2*psi) \
                + 2 * cos(iota) * cos (phi) * cos (2*psi) ) ]

    return a


# Define the time dependent basis functions as shown in the equation (19)
def Basis (theta, phi, theta_p, phi_p, t):
    A = []
    F = AntennaPattern(theta, phi, theta_p, phi_p, t)
    A += [ F[0] * omega(t)**(-1./3) * sin (2 * Phi(t)) ]
    A += [ F[0] * omega(t)**(-1./3) * cos (2 * Phi(t)) ]
    A += [ F[1] * omega(t)**(-1./3) * sin (2 * Phi(t)) ]
    A += [ F[1] * omega(t)**(-1./3) * cos (2 * Phi(t)) ]

    return A

# Define the pulsar term as in the eq (17)
#def pulsar (t, zeta, iota, Phi0 psi, theta, phi, omega0, L):
