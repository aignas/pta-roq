
#include <cmath>
#include <vector>
#include <Vec3.hh>

// Antena Pattern function
std::vector<double> antennaPattern (double theta, double phi, cav::Vec3 u_p) {
    // Initialise unit vector according to the given parameters
    UnitVectors v (theta, phi)

    double d1, d2, d3
    std::vector<double> F

    d1 = cav.Vec3.dot(v.m(),u_p)
    d2 = cav.Vec3.dot(v.n(),u_p)
    d3 = cav.Vec3.dot(v.Omega(),u_p)

    F[0] = 0.5 * ( d1 + d2 ) * ( d1 - d2 ) / ( 1 + d3 )
    F[1] = ( d1 * d2 ) / ( 1 + d3 )

    return F
}

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
    cdef int i

    # Add contributions to the signal from all GW sources
    for i in range(0,sources.shape[0],8):
        Signal += individualSource(t, sources[i:i+8], L, u_p)

    # Add noise to the signal
    Signal += noise (t, variance)

    return Signal

def genSchedule (np.ndarray schedule, double t_final, double dt_min, double dt_max):
    cdef double t
    cdef np.ndarray u_p
    cdef np.ndarray dates_out = np.array([]), index_out = np.array([])
    cdef int N = schedule.shape[0]
    collectData = True
    cdef int i

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
    cdef int i

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
    cdef int k

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
    for k in range(N):
        sum_member *= - (f * tau)**2 / ((2*k + 1) * (2*k + 2)) * (2*k - 1 - gamma) \
                      / (2*k + 1 - gamma)
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
    cdef int k

    if a==b:
        # Here I use a similar technique to the one explained in
        # covarianceMatrixMemberGWB
        sum_member = 1 / (1 - gamma)
        sum = sum_member
        for k in range(N):
            sum_member *= - (f * tau)**2 / ((2*k + 1) * (2*k + 2)) * (2*k + 1 - gamma) \
                          / (2*k + 3 - gamma)
            sum += sum_member

        r =  A**2 / (f**(gamma - 1)) \
             * (np.gamma(1 - gamma) * sin(np.pi * gamma /2) * (f * tau)**( gamma -1) - sum)
    
    # Return the computed value
    return r
