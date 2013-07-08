#include <vector>
#include <valarray>
#include <cmath>

#include "linalg.hh"
#include "pulsar.hh"
#include "signal-model.hh"
#include "signal-sampling.hh"

void generateSample (std::vector<double>& out, std::vector<Pulsar> &pulsars,
                     std::vector<unsigned short>& indices, std::vector<double>& Times,
                     std::vector<std::vector<double> > &sources) { 
    Pulsar pulsar;
    const unsigned int N = indices.size();

    if (out.size() != N) {
        out.resize(N);
    }

    // Start collecting the data
    for (unsigned int i = 0; i < N; i++) {
        out[i] = residual(Times[i], indices[i], sources, pulsar);
    }
}

void genCovarianceMatrix (std::vector<double> & matrix, std::vector<Pulsar> & pulsars,
                          std::vector<unsigned short> & indices, std::vector<double> & Times,
                          const bool WhiteNoise, const bool RedNoise, const bool PowerLaw, const bool GWB) {
    unsigned int N = indices.size();

    // Initiate a zero matrix
    matrix.resize(N*N);

    // Set the low frequency cutoff by dividing frequency error by the time of the
    // measurement
    double f_L = 1e-5 / Times.back();

    // Initiate the sum truncation in GWB and Power Law Noise terms
    unsigned int NTrunk = 10000;

    if (WhiteNoise) {
        unsigned int pidx = 0;
        double amplitude = pow(pulsars[pidx].getWhiteNoise(), 2);
        for (unsigned int i = 0; i < N; i++) {
            if ( pidx != indices[i]) {
                pidx = indices[i];
                amplitude = pow(pulsars[pidx].getWhiteNoise(), 2);
            }

            // Fill only the diagonal entries
            matrix[ i * (N+1) ] = amplitude;
        }
    }

    if (RedNoise) {
    }

    if (PowerLaw) {
    }

    if (GWB) {
        // Assume that the noise comes from SMBHBs only
        double gamma = 7./3;

        // Let the GWB amplitude be small
        double A_GWB = 1e-15, tau, C;
        std::vector<double> uPV1 (3), uPV2 (3);

        unsigned int a = 0, b = 0;
        for (unsigned int i = 0; i < N; i++) {
            for (unsigned int j = 0; j < N; j++) {
                a = indices[i];
                b = indices[j];

                // Calculate the difference in times
                tau = 2 * _M_PI * ( Times[i] - Times[j] );

                // From the pulsars we need only the unit vectors and the noise
                // magnitude, we could probably do it here, as it would make the code
                // look better
                uPV1 = pulsars[a].getUnitVector();
                uPV2 = pulsars[b].getUnitVector();
                C = 1./2 * (1 - dotProduct( uPV1, uPV2 ));

                // Fill only the diagonal entries
                matrix[ i*j ] = covarianceMatrixMemberGWB (i, j, a, b, A_GWB, f_L, gamma, tau, NTrunk, C);

            }
        }
    }

    // Invert the matrix, as we do not need the other bit
    inverse(matrix);
}

// Calculate red noise Lorentzian terms (i.e. red noise)
double covarianceMatrixMemberLor (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double N, double f, double tau){
    
    // r is the return value and N should be the noise amplitude
    double r = 0;
    
    if (a==b) {
        r = N * N * exp(-f*tau);
    }

    return r;
}

// Calculate the power law spectral noise
double covarianceMatrixMemberPowLaw (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double A, double f, double tau, double gamma, unsigned int N) {
    
    double r = 0;

    if (a==b) {
        // Here I use a similar technique to the one explained in
        // covarianceMatrixMemberGWB
        double sum_member = 1 / (1 - gamma);
        double sum = sum_member;
        for (unsigned int k = 0; k < N; k++ ) {
            sum_member *= - pow((f * tau), 2) / ((2*k + 1) * (2*k + 2)) * (2*k + 1 - gamma) \
                          / (2*k + 3 - gamma);
            sum += sum_member;
        }

        r =  A * A / pow(f, gamma - 1) \
             * (tgamma(1 - gamma) * sin(_M_PI * gamma /2) * pow(f * tau, gamma - 1) - sum);
    }
    
    // Return the computed value
    return r;
}


// Calculate the GWB terms in the covariance matrix members.
double covarianceMatrixMemberGWB (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double A, double f, double gamma, double tau, unsigned int N, double C) {

    double alpha = 3/2 * C * log(C) - C/4 + 1/2;
    // a simple delta function implementation
    if (a == b) {
        alpha += 1/2;
    }

    // Here I use slightly more memory by storring each_member before summing it, but
    // this way I do not have to calculate horrible factorials and it should speed things
    // up a bit
    // This function calculates N terms and then truncates the series
    double sum_member = - 1 / (1 + gamma);
    double sum = sum_member;
    for (unsigned int k = 0; k < N; k++) {
        sum_member *= - pow(f * tau, 2) / ((2*k + 1) * (2*k + 2)) * (2*k - 1 - gamma) \
                      / (2*k + 1 - gamma);
        sum += sum_member;
    }

    return A * A * alpha / (pow(2 * _M_PI, 2) * pow(f, 1 + gamma)) \
           / (tgamma(-1 - gamma) * sin(-_M_PI * gamma / 2) * pow(f * tau, 1 + gamma) - sum);
}
