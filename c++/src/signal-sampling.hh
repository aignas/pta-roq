/**
 * This is a header file for the signal sampler
 */

#include <vector>
#include <valarray>
#include "pulsar.hh"

#ifndef _SIGNAL_SAMPLING_HXX_
#define _SIGNAL_SAMPLING_HXX_
// Define PI to a high precision
const double _M_PI = 3.141592653589793238462643383279502884;
#endif

/**
 * Generate a signal, with optionally specifying whether to include noise or not
 *
 * @param out The output array
 * @param schedule The timing schedule, which should be use for data generation.
 * @param pulsars A list of pulsars
 * @param indices The index of pulsars in the schedule
 * @param Times The time series
 * @param sources A list of sources
 */
void generateSample (std::vector<double> & out, std::vector<Pulsar> & pulsars, 
                     std::vector<unsigned short> & indices, std::vector<double> & Times,
                     std::vector<std::vector<double> > & sources);

/**
 * Generate a covariance matrix, which depends on the sampling parameters (i.e. the
 *  schedule
 *
 * @param matrix_out An array, which will store the matrix
 * @param pulsars A grid of pulsars
 * @param WhiteNoise If true, then white noise contributions will be added.
 * @param RedNoise If true, then red noise contributions will be added.
 * @param PowerLaw If true, then power law noise contributions will be added.
 * @param GWB If true, then GWB contributions will be added.
 */
void genCovarianceMatrix (std::vector<double> &matrix, std::vector<Pulsar> & pulsars,
                          std::vector<unsigned short> & indices, std::vector<double> & Times,
                          const bool WhiteNoise, const bool RedNoise, const bool PowerLaw, const bool GWB); 

/**
 * Calculate the GWB terms in the covariance matrix members.
 *
 * @param i measurement index
 * @param j measurement index
 * @param a pulsar index
 * @param b pulsar index
 * @param A GWB amplitude
 * @param f cutoff frequency
 * @param gamma Some statistical coefficient describing the population of the sources of
 *    GWB
 * @param tau characteristic time
 * @param N The number of terms in the sum
 * @param C some numerical factor, which is needed in the formula
 */
double covarianceMatrixMemberGWB (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double A, double f, double gamma, double tau, unsigned int N, double C);

/**
 * Calculate the power law spectral noise
 *
 * @param i measurement index
 * @param j measurement index
 * @param a pulsar index
 * @param b pulsar index
 * @param A Power Law noise amplitude
 * @param f cutoff frequency
 * @param tau characteristic time
 * @param gamma Some statistical coefficient describing the power law noise
 * @param N The number of terms in the sum
 */
double covarianceMatrixMemberPowLaw (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double A, double f, double tau, double gamma, unsigned int N);

/**
 * Calculate red noise Lorentzian terms (i.e. red noise)
 *
 * @param i measurement index
 * @param j measurement index
 * @param a pulsar index
 * @param b pulsar index
 * @param N Red noise amplitude
 * @param f characteristic frequency
 * @param tau characteristic time
 */
double covarianceMatrixMemberLor (unsigned int i, unsigned int j, unsigned int a, unsigned int b, 
        double N, double f, double tau);
