#include <vector>
#include <pulsar.hh>
#include <signal-vectors.hh>
#include <signal-model.hh>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_rng.h>

#include <time.h>

random (double LO, double HI) {
    // FIXME: I have an example. Look it up
    gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
    for ( int i = 0; i < 5; i++ )
        std::cout << gsl_rng_uniform( rng ) << ", ";
    std::cout << "\n";
    gsl_rng_free( rng );

    return LO + (double)rand()/((double)RAND_MAX/(HI-LO));
}

dvec antennaPattern (dvec intrinsic, dvec u_p) {
    // Using at method makes checks whether we are out of bounds. It is just better to
    // be sure.
    UnitVectors v (intrinsic.at(0), intrinsic.at(1));

    // Initialise some parameters
    double d1 = 0, 
           d2 = 0, 
           d3 = 0;
    dvec F (2);

    // Calculate three dot products
    d1 = dotProduct(v.m(), u_p);
    d2 = dotProduct(v.n(), u_p);
    d3 = dotProduct(v.Omega(), u_p);

    F[0] = 0.5 * ( d1 + d2 ) * ( d1 - d2 ) / ( 1 + d3 );
    F[1] = ( d1 * d2 ) / ( 1 + d3 );

    return F;
}

double omegaReduced(double t, double omega0, double M) {
    return pow(pow(omega0, -8./3) - 256./5 * pow(M, 5./3) * t, 1./8);
}

double Phi(double, double, double) {
    return 1./(32 * pow(M, 5./3)) * ( pow(omega0, -5./3) - pow(omegaReduced(t, omega0, M), 5));
}

std::vector<double> gravWaveContrib (double t, dvec i, double omega) {
    // Assign the intrinsic parameters to the names
    double M    = i.at(0),
           D    = i.at(1),
           iota = i.at(2),
           Phi0 = i.at(3),
           psi  = i.at(4);

    // Create some temp vars
    double a, b, zeta = M/D;
    dvec gw (2);

    a = sin( 2 * (Phi(t, omega0, M) - Phi0) ) * (1 + pow(cos(iota), 2));
    b = cos( 2 * (Phi(t, omega0, M) - Phi0) ) * cos(iota);

    gw[0] = zeta * omegaReduced(t, omega0, M) * ( - a * cos(2*psi) - 2 * sin(2*psi) );
    gw[1] = zeta * omegaReduced(t, omega0, M) * ( - a * sin(2*psi) - 2 * cos(2*psi) );

    return gw;
}

dvec amplitude (dvec i) {
    // Assign the intrinsic parameters to the names
    double zeta = i.at(0) / intrinsic.at(1),
           iota = i.at(2),
           Phi0 = i.at(3),
           psi  = i.at(4);

    // Assign some temp variables:
    double c = cos(iota);
    double b = 1 + pow(c, 2);

    dvec a (4);

    a[0] =  zeta * ( b * cos (phi) * cos (2*psi) + 2 * c * sin (phi) * sin (2*psi) );
    a[1] = -zeta * ( b * sin (phi) * cos (2*psi) + 2 * c * cos (phi) * sin (2*psi) );
    a[2] =  zeta * ( b * cos (phi) * sin (2*psi) + 2 * c * sin (phi) * cos (2*psi) );
    a[3] = -zeta * ( b * sin (phi) * sin (2*psi) + 2 * c * cos (phi) * cos (2*psi) );

    return a;
}

dvec basis (double t, dvec extrinsic, dvec u_p) {
    // Assign the extrinsic parameters to the names
    double theta  = extrinsic.at(0),
           phi    = extrinsic.at(1),
           omega0 = extrinsic.at(2);

    dvec A (4), 
         F = antennaPattern (theta, phi, u_p);

    A[0] = F.at(0) * omegaReduced(t, omega0, M) * sin (2 * Phi(t, omega0, M));
    A[1] = F.at(0) * omegaReduced(t, omega0, M) * cos (2 * Phi(t, omega0, M));
    A[2] = F.at(1) * omegaReduced(t, omega0, M) * sin (2 * Phi(t, omega0, M));
    A[3] = F.at(1) * omegaReduced(t, omega0, M) * cos (2 * Phi(t, omega0, M));

    return A;
}

double pulsarTerm (double t, dvec i, dvec e, double L, dvec pUV) {
    // unpack the parameters
    double M      = i.at(0),
           D      = i.at(1),
           iota   = i.at(2),
           Phi0   = i.at(3),
           psi    = i.at(4),
           theta  = e.at(0),
           phi    = e.at(1),
           omega0 = e.at(2);

    // Initialize the unit vector data struct
    UnitVectors v (intrinsic.at(0), intrinsic.at(1));

    // Calculate the time at the pulsar
    double tp = t - L * (1 + dotProduct(v.Omega(), u_p));

    // Calculate
    dvec F = antennaPattern(theta, phi, pUV),
         s = gWContribution(tp, intrinsic, omega0);

    return dotProduct(F,s);
}

double noise (double t, double whiteNoise, dvec redNoise, dvec powerLawNoise) {
    // FIXME: implement other types of noise as well
    double white = random() * pow(whiteNoise, 0.5),
           red = 0,
           powerLaw = 0;

    return white + red + powerLaw;
}

double individualSource (double t, dvec intrinsic, dvec extrinsic, double L, dvec pUV) {

    dvec a = amplitude(intrinsic),
         A = basis (t, extrinsic, pUV);

    double p = (t, intrinsic, extrinsic, L, pUV);

    return dotProduct(a,A) + p;
}

// FIXME: Here I need a structure, so that I would have only a single pulsar properties
double residual (double t, unsigned int N, std::vector<dvec> sources, dvec pulsarProperties) {
    // Separate the pulsar properties into the unit Vector, distance and noise
    // properties

    // Add contributions to the signal from all GW sources
    double Signal = 0;
    for (unsigned int i = 0; i < sources.size(); i++) {
        Signal += individualSource(t, sources[i], L, u_p)
    }

    // Add noise to the signal
    Signal += noise (t, white, red, pwlaw);

    return Signal;
}
