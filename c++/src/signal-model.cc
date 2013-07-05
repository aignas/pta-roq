#include <vector>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_rng.h>

#include <time.h>

#include "pulsar.hh"
#include "signal-vectors.hh"
#include "signal-model.hh"

// load the linear algebra helpers
#include "linalg.hh"

double random_gsl (double LO, double HI) {
    // FIXME: I have an example. Look it up
    gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_uniform( rng );

    // Free the generator
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

double omegaReduced (double omega0) {
    // Make a simplification, which is sensible like in [Ellis et al. 2010]_
    //return pow(pow(omega0, -8./3) - 256./5 * pow(M, 5./3) * t, 1./8);
    return pow(omega0, -1./3);
}

double Phi(double t, double omega0) {
    // Make a simplification here as well
    //return 1./(32 * pow(M, 5./3)) * ( pow(omega0, -5./3) - pow(omegaReduced(t, omega0, M), 5));
    return pow(omegaReduced (omega0), -3) * t;
}

std::vector<double> gravWaveContrib (double t, dvec i, double omega0) {
    // Assign the intrinsic parameters to the names
    double M    = i.at(0),
           D    = i.at(1),
           iota = i.at(2),
           Phi0 = i.at(3),
           psi  = i.at(4);

    // Create some temp vars
    double a, b, zeta = M/D;
    dvec gw (2);

    a = sin( 2 * (Phi(t, omega0) - Phi0) ) * (1 + pow(cos(iota), 2));
    b = cos( 2 * (Phi(t, omega0) - Phi0) ) * cos(iota);

    gw[0] = zeta * omegaReduced(omega0) * ( - a * cos(2*psi) - 2 * sin(2*psi) );
    gw[1] = zeta * omegaReduced(omega0) * ( - a * sin(2*psi) - 2 * cos(2*psi) );

    return gw;
}

dvec amplitude (dvec i) {
    // Assign the intrinsic parameters to the names
    double zeta = i.at(0) / i.at(1),
           iota = i.at(2),
           Phi0 = i.at(3),
           psi  = i.at(4);

    // Assign some temp variables:
    double c = cos(iota);
    double b = 1 + pow(c, 2);

    dvec a (4);

    a[0] =  zeta * ( b * cos (Phi0) * cos (2*psi) + 2 * c * sin (Phi0) * sin (2*psi) );
    a[1] = -zeta * ( b * sin (Phi0) * cos (2*psi) + 2 * c * cos (Phi0) * sin (2*psi) );
    a[2] =  zeta * ( b * cos (Phi0) * sin (2*psi) + 2 * c * sin (Phi0) * cos (2*psi) );
    a[3] = -zeta * ( b * sin (Phi0) * sin (2*psi) + 2 * c * cos (Phi0) * cos (2*psi) );

    return a;
}

dvec basis (double t, dvec extrinsic, dvec u_p) {
    // Assign the extrinsic parameters to the names
    double theta  = extrinsic.at(0),
           phi    = extrinsic.at(1),
           omega0 = extrinsic.at(2);

    dvec A, F = antennaPattern (extrinsic, u_p);

    A.push_back(F.at(0) * omegaReduced(omega0) * sin (2 * Phi(t, omega0)));
    A.push_back(F.at(0) * omegaReduced(omega0) * cos (2 * Phi(t, omega0)));
    A.push_back(F.at(1) * omegaReduced(omega0) * sin (2 * Phi(t, omega0)));
    A.push_back(F.at(1) * omegaReduced(omega0) * cos (2 * Phi(t, omega0)));

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
    UnitVectors v (theta, phi);

    // Calculate the time at the pulsar
    double tp = t - L * (1 + dotProduct(v.Omega(), pUV));

    // Calculate
    dvec F = antennaPattern(e, pUV),
         s = gravWaveContrib(tp, i, omega0);

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

    double p = pulsarTerm (t, intrinsic, extrinsic, L, pUV);

    return dotProduct(a,A) + p;
}

// FIXME: Here I need a structure, so that I would have only a single pulsar properties
double residual (double t, unsigned int N, std::vector<dvec> sources, Pulsar pulsarProperties) {
    // Separate the pulsar properties into the unit Vector, distance and noise
    // properties

    // Add contributions to the signal from all GW sources
    double Signal = 0, 
           L = pulsarProperties.getDistance(),
           white = pulsarProperties.getWhiteNoise();
    dvec pUV = pulsarProperties.getUnitVector(),
         red = pulsarProperties.getRedNoise(),
         pwlaw = pulsarProperties.getPowerLawNoise();

    for (unsigned int i = 0; i < sources.size(); i++) {
        // FIXME Think whether I need an additional data structure here
        Signal += individualSource(t, sources[i], L, pUV);
    }

    // Add noise to the signal
    Signal += noise (t, white, red, pwlaw);

    return Signal;
}
