#include <vector>
#include <cmath>
#include <iostream>

#include "vectors.hh"
#include "model.hh"
#include "../pulsar.hh"
#include "../random-helper.hh"

// load the linear algebra helpers
#include "../linalg.hh"

void antennaPattern (dvec& intrinsic, dvec& u_p, dvec& out) {
    // Using at method makes checks whether we are out of bounds. It is just better to
    // be sure.
    UnitVectors v (intrinsic.at(0), intrinsic.at(1));
    // Initialise some parameters
    dvec vm = v.m(),
         vn = v.n(),
         vOmega = v.Omega();

    // Calculate three dot products
    const double d1 = dotProduct(vm, u_p),
                 d2 = dotProduct(vn, u_p),
                 d3 = dotProduct(vOmega, u_p);

    out.resize(2);
    out[0] = 0.5 * ( d1 + d2 ) * ( d1 - d2 ) / ( 1 + d3 );
    out[1] = ( d1 * d2 ) / ( 1 + d3 );
}

inline double omegaReduced (const double omega0) {
    // Make a simplification, which is sensible like in [Ellis et al. 2010]_
    //return pow(pow(omega0, -8./3) - 256./5 * pow(M, 5./3) * t, 1./8);
    return pow(omega0, -1./3);
}

inline double PhiDash(const double t, const double omega0) {
    // Make a simplification here as well
    //return 1./(32 * pow(M, 5./3)) * ( pow(omega0, -5./3) - pow(omegaReduced(t, omega0, M), 5));
    return pow(omegaReduced (omega0), -3) * t;
}

void amplitude (dvec& e, dvec& out) {
    // Assign the intrinsic parameters to the names
    const double zeta = e.at(0),
                 iota = e.at(1),
                 phi0 = e.at(2),
                 psi  = e.at(3),
                 // Assign some temp variables:
                 c = cos(iota),
                 b = 1 + pow(c, 2);

    out.resize(4);
    out[0] =  zeta * ( b * cos (2*phi0) * cos (2*psi) + 2 * c * sin (2*phi0) * sin (2*psi) );
    out[1] = -zeta * ( b * sin (2*phi0) * cos (2*psi) + 2 * c * cos (2*phi0) * sin (2*psi) );
    out[2] =  zeta * ( b * cos (2*phi0) * sin (2*psi) + 2 * c * sin (2*phi0) * cos (2*psi) );
    out[3] = -zeta * ( b * sin (2*phi0) * sin (2*psi) + 2 * c * cos (2*phi0) * cos (2*psi) );
}

void basis (const double t, dvec& intrinsic, dvec& u_p, dvec& out, double delta) {
    // Assign the intrinsic parameters to the names
    const double omega0 = intrinsic.at(2),
                 // The following expressions can only be used in the approximation,
                 // that the frequency of the signal from the SMBHB doesn't change.
                 // Otherwise it is not possible.
                 a = sin(2 * PhiDash(t, omega0)),
                 b = sin(2 * PhiDash(delta, omega0)),
                 c = cos(2 * PhiDash(t, omega0)),
                 d = cos(2 * PhiDash(delta, omega0));

    dvec F;
    antennaPattern (intrinsic, u_p, F);

    // This fixes some numerical issues with sin(i+j) formula if i >> j
    out.resize(4);
    out[0] = F.at(0) * omegaReduced(omega0) * (a * d - c * b);
    out[1] = F.at(0) * omegaReduced(omega0) * (c * d + a * b); 
    out[2] = F.at(1) * omegaReduced(omega0) * (a * d - c * b); 
    out[3] = F.at(1) * omegaReduced(omega0) * (c * d + a * b); 
}

double noise (const double t, const double whiteNoise, dvec& redNoise, dvec& powerLawNoise) {
    // FIXME: implement other types of noise as well
    double white = random_uniform(0, pow(whiteNoise, 0.5)),
           red = 0,
           powerLaw = 0;

    random_uniform_free();

    return white + red + powerLaw;
}

double individualSource (const double t, dvec& params, const double L, dvec& pUV, bool includePulsarTerm) {
    dvec extrinsic (4, 0),
         intrinsic (3, 0),
         a, A, A_p;

    // Translate params
    extrinsic[0] = params.at(0);
    extrinsic[1] = params.at(1);
    extrinsic[2] = params.at(2);
    extrinsic[3] = params.at(3);
    intrinsic[0] = params.at(4);
    intrinsic[1] = params.at(5);
    intrinsic[2] = params.at(6);

    if (fabs(intrinsic[0]) < 1e-9 or fabs(intrinsic[0] - _M_PI) < 1e-9) {
        throw "Error: The azimuthal angle is close to the coordinate singularity";
    } else if (fabs(intrinsic[1]) > _M_PI) {
        throw "Error: The polar angle is out of range (-pi;pi)";
    }
    // TODO other things to check
    // omega out of detection range
    // mass out of detection range
    // iota, psi out of range

    // Calculate the time at the pulsar
    UnitVectors v (intrinsic.at(0), intrinsic.at(1));
    dvec vOmega = v.Omega();
    double delta = L * (1 + dotProduct(vOmega, pUV));

    amplitude(extrinsic, a);
    basis(t, intrinsic, pUV, A);
    basis(t, intrinsic, pUV, A_p, delta);

    // Do not use the pulsarTerm function, as it is giving nonsenses because tp >> 1,
    // which means that tp + something will be not precise when we take a sine of it.
    // This is now implemented with the additional (optional) parameter in the basis
    // function.
    if (includePulsarTerm) {
        axpyProduct(-1,A_p,A);
    }

    return dotProduct(a,A);
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

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Functions which should not be used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// NOTE This function doesn't give correct results at t >> 1 because of the phidash
// approximation. This function should not be used!!!!
void gravWaveContrib (const double t, dvec& e, const double omega0, dvec& out) {
    // Assign the intrinsic parameters to the names
    const double zeta = e.at(0),
                 iota = e.at(1),
                 phi0 = e.at(2),
                 psi  = e.at(3);

    // Create some temp vars
    const double a = - sin( 2 * (PhiDash(t, omega0) - phi0) ) * (1 + pow(cos(iota), 2)),
                 b = - 2 * cos( 2 * (PhiDash(t, omega0) - phi0)) * cos(iota);

    out.resize(2);
    out[0] = zeta * omegaReduced(omega0) * ( a * cos(2*psi) +  b * sin(2*psi) );
    out[1] = zeta * omegaReduced(omega0) * ( a * sin(2*psi) -  b * cos(2*psi) );
}

// Since this uses the gravWaveContrib, this shouldn't be used either!!!!!
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
    dvec F, s;
    antennaPattern(i, pUV, F);
    gravWaveContrib(tp, e, omega0, s);

    return dotProduct(F,s);
}
