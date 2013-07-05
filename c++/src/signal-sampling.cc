#include <vector>
#include <cmath>

#include "pulsar.hh"
#include "signal-sampling.hh"
#include "signal-model.hh"
#include "linalg.hh"

// Define PI to a high precision
const static double M_PI = 3.141592653589793238462643383279502884;

std::vector<double> generateSample (PulsarGrid pulsars, std::vector<std::vector<double> > sources) {
    double t, L;
    std::vector<double> u_p, noise, res, Times;
    unsigned int N = pulsars.getNumber();


    // Start collecting the data
    for (unsigned int i = 0; i < N; i++) {
        Times = pulsars.getSchedule(i);
        for (unsigned int j = 0; j < Times.size(); j++) {
            L = pulsars.getDistance(i);
            u_p = pulsars.getUnitVector(i);
            t = Times[j];

            Pulsar pulsarProps = pulsars.getPulsarProperties(i);

            res.push_back(residual(t, N, sources, pulsarProps));
        }
    }

    return res;

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
        double A, double f, double tau, double gamma, unsigned int N){
    
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
             * (tgamma(1 - gamma) * sin(M_PI * gamma /2) * pow(f * tau, gamma - 1) - sum);
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

    return A * A * alpha / (pow(2 * M_PI, 2) * pow(f, 1 + gamma)) \
           / (tgamma(-1 - gamma) * sin(-M_PI * gamma / 2) * pow(f * tau, 1 + gamma) - sum);
}

std::valarray<double>* covarianceMatrix (std::vector<Pulsar> pulsars, 
        bool WhiteNoise, bool RedNoise, bool PowerLaw, bool GWB) {

    unsigned int N = 0, Nt;
    std::vector<double> timestamps;
    for (unsigned int i = 0; i < pulsars.size(); i++) {
        timestamps.push_back(pulsars[i].getSchedule().size());
        N += timestamps.at(i);
    }

    // Initiate a zero matrix
    std::valarray<double> matrix (N*N);

    // Set the low frequency cutoff by dividing frequency error by the time of the
    // measurement
    double f_L = 1e-5 / pulsars[0].getSchedule().back();

    // Initiate the sum truncation in GWB and Power Law Noise terms
    unsigned int NTrunk = 10000;

    if WhiteNoise {
        // Pulsar id
        unsigned int pidx = 0, j=0;
        double amplitude = pow(pulsars[pidx].getWhiteNoise(), 2);
        for (unsigned int i = 0; i < N; i++) {
            // Keep track of the pulsar index with comparing the secondary counter with
            // the number of timestamps for that pulsar.
            if (timestamps[pidx] == j) {
                pidx++; j = 0;
                amplitude = pow(pulsars[pidx].getWhiteNoise(), 2);
            } else {
                j++;
            }

            // Fill only the diagonal entries
            matrix( i * (N+1) ) = amplitude;
        }
    }

    if RedNoise {
    }

    if PowerLawNoise {
    }

    if GWB {
        // Assume that the noise comes from SMBHBs only
        double gamma = 7./3;

        // Let the GWB amplitude be small
        double A_GWB = 1e-15,
               a, b, tau, C;

        unsigned int a = 0, b = 0, k = 0, l = 0;
        double amplitude = pow(pulsars[pidx].getWhiteNoise(), 2);
        for (unsigned int i = 0; i < N; i++) {
            // Keep track of the pulsar index with comparing the secondary counter with
            // the number of timestamps for that pulsar. This is very similar to what I
            // did in the white noise case
            if (timestamps[a] == k) {
                a++; k = 0;
            } else {
                k++;
            }

            for (unsigned int j = 0; j < N; j++) {
                if (timestamps[pidx] == j) {
                    b++; l = 0;
                } else {
                    l++;
                }
                // Calculate the difference in times
                tau = 2 * M_PI * (pulsars[a].getSchedule()[k] - pulsars[b].getSchedule()[l]);

                // From the pulsars we need only the unit vectors and the noise
                // magnitude, we could probably do it here, as it would make the code
                // look better
                C = 1./2 * (1 - dotProduct( pulsars[a].getUnitVector(), pulsars[b].getUnitVector()))

                // Fill only the diagonal entries
                matrix( i*j ) = covarianceMatrixMemberGWB (i, j, a, b, A_GWB, f_L, gamma, tau, NTrunk, C);

            }
        }
    }

    // Return by reference, otherwise, we will have lots of copying
    return &matrix;
}

