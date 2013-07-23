#include <vector>
#include <cmath>
#include <iostream>

#include "vectors.hh"
#include "model.hh"
#include "../pulsar.hh"
#include "../random-helper.hh"

// load the linear algebra helpers
#include "linalg.hh"

dvec antennaPattern (dvec& intrinsic, dvec& u_p) {
    // Using at method makes checks whether we are out of bounds. It is just better to
    // be sure.
    UnitVectors v (intrinsic.at(0), intrinsic.at(1));
    dvec vm = v.m(),
         vn = v.n(),
         vOmega = v.Omega();

    // Initialise some parameters
    dvec F (2);

    // Calculate three dot products
    double d1 = dotProduct(vm, u_p),
           d2 = dotProduct(vn, u_p),
           d3 = dotProduct(vOmega, u_p);

    F[0] = 0.5 * ( d1 + d2 ) * ( d1 - d2 ) / ( 1 + d3 );
    F[1] = ( d1 * d2 ) / ( 1 + d3 );

    return F;
}

inline double omegaReduced (const double omega0) {
    // Make a simplification, which is sensible like in [Ellis et al. 2010]_
    //return pow(pow(omega0, -8./3) - 256./5 * pow(M, 5./3) * t, 1./8);
    return pow(omega0, -1./3);
}

inline double Phi(const double t, const double omega0) {
    // Make a simplification here as well
    //return 1./(32 * pow(M, 5./3)) * ( pow(omega0, -5./3) - pow(omegaReduced(t, omega0, M), 5));
    return pow(omegaReduced (omega0), -3) * t;
}

std::vector<double> gravWaveContrib (const double t, dvec& e, const double omega0) {
    // Assign the intrinsic parameters to the names
    double M    = e.at(0),
           D    = e.at(1),
           iota = e.at(2),
           Phi0 = e.at(3),
           psi  = e.at(4);

    // Create some temp vars
    double zeta = M/D;
    dvec gw (2);

    double a = - sin( 2 * (Phi(t, omega0) - Phi0) ) * (1 + pow(cos(iota), 2)),
           b = - 2 * cos( 2 * (Phi(t, omega0) - Phi0) ) * cos(iota);

    gw[0] = zeta * omegaReduced(omega0) * ( a * cos(2*psi) +  b * sin(2*psi) );
    gw[1] = zeta * omegaReduced(omega0) * ( a * sin(2*psi) -  b * cos(2*psi) );

    return gw;
}

dvec amplitude (dvec& e) {
    // Assign the intrinsic parameters to the names
    double zeta = e.at(0) / e.at(1),
           iota = e.at(2),
           Phi0 = e.at(3),
           psi  = e.at(4);

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

dvec basis (const double t, dvec& intrinsic, dvec& u_p) {
    // Assign the intrinsic parameters to the names
    double omega0 = intrinsic.at(2);

    dvec A, F = antennaPattern (intrinsic, u_p);

    A.push_back(F.at(0) * omegaReduced(omega0) * sin (2 * Phi(t, omega0)));
    A.push_back(F.at(0) * omegaReduced(omega0) * cos (2 * Phi(t, omega0)));
    A.push_back(F.at(1) * omegaReduced(omega0) * sin (2 * Phi(t, omega0)));
    A.push_back(F.at(1) * omegaReduced(omega0) * cos (2 * Phi(t, omega0)));

    return A;
}

double pulsarTerm (const double t, dvec& e, dvec& i, const double L, dvec& pUV) {
    // unpack the parameters
    double theta  = i.at(0),
           phi    = i.at(1),
           omega0 = i.at(2);

    // Initialize the unit vector data struct
    UnitVectors v (theta, phi);
    dvec vOmega = v.Omega();

    // Calculate the time at the pulsar
    double tp = t - L * (1 + dotProduct(vOmega, pUV));

    // Calculate
    dvec F = antennaPattern(i, pUV),
         s = gravWaveContrib(tp, e, omega0);

    return dotProduct(F,s);
}

double noise (const double t, const double whiteNoise, dvec& redNoise, dvec& powerLawNoise) {
    // FIXME: implement other types of noise as well
    double white = random_uniform(0, pow(whiteNoise, 0.5)),
           red = 0,
           powerLaw = 0;

    random_uniform_free();

    return white + red + powerLaw;
}

double individualSource (const double t, dvec& params, const double L, dvec& pUV) {
    dvec extrinsic (5, 0),
         intrinsic (3, 0);

    // Translate params
    extrinsic[0] = params.at(0);
    extrinsic[1] = params.at(1);
    extrinsic[2] = params.at(2);
    extrinsic[3] = params.at(3);
    extrinsic[4] = params.at(4);
    intrinsic[0] = params.at(5);
    intrinsic[1] = params.at(6);
    intrinsic[2] = params.at(7);

    if (fabs(extrinsic[0]) <  1e-9) {
        throw "Error: The chirp mass of the source is too low";
    } else if (fabs(extrinsic[1]) < 1e-9) {
        throw "Error: The distance to the source is too low";
    } else if (fabs(intrinsic[0]) < 1e-9 or fabs(intrinsic[0] - _M_PI) < 1e-9) {
        throw "Error: The azimuthal angle is close to the coordinate singularity";
    } else if (fabs(intrinsic[1]) > _M_PI) {
        throw "Error: The polar angle is out of range (-pi;pi)";
    }
    // TODO other things to check
    // omega out of detection range
    // mass out of detection range
    // iota, psi out of range

    dvec a = amplitude(extrinsic),
         A = basis (t, intrinsic, pUV);

    // FIXME Do not use the pulsar term for the time being. Think of a switch to implement it
    //double p = pulsarTerm (t, extrinsic, intrinsic, L, pUV);
    double p = 0;

    return dotProduct(a,A) + p;
}

// FIXME: Here I need a structure, so that I would have only a single pulsar properties
double residual (const double t, const unsigned int N, std::vector<dvec>& sources, Pulsar& pulsarProperties, bool noiseInclude) {
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
        try {
            Signal += individualSource(t, sources[i], L, pUV);
        } catch (const char* mes) {
            std::cerr << mes << std::endl;
        }
    }

    // Add noise to the signal
    if (noiseInclude) {
        Signal += noise (t, white, red, pwlaw);
    }

    return Signal;
}
