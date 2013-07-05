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
void generateSample (std::vector<double>& out, std::vector<Pulsar> &pulsars, 
                     std::vector<unsigned short>& indices, std::vector<double>& Times,
                     std::vector<std::vector<double> > &sources);

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
void covarianceMatrix (std::valarray<double> &matrix_out, std::vector<Pulsar> pulsars, 
        bool WhiteNoise = true, bool RedNoise = true, bool PowerLaw = true, bool GWB = true);
