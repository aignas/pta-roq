/**
 * This is a header file for the signal sampler
 */

#include <vector>
#include "pulsar.hh"

/**
 * Generate a signal, with optionally specifying whether to include noise or not
 *
 * @param schedule The timing schedule, which should be use for data generation.
 * @param pulsars A list of pulsars
 * @param sources A list of sources
 */
std::vector<double> generateSample (PulsarGrid pulsars, std::vector<std::vector<double>> sources);

/**
 * Generate a covariance matrix, which depends on the sampling parameters (i.e. the
 *  schedule
 */
std::valarray<double>* covarianceMatrix (std::vector<Pulsar> pulsars, 
        bool whiteNoise = true, bool RedNoise = false, bool PowerLaw = false, bool GWB = false) {
