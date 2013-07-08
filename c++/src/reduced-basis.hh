
#include <vector>
#include <valarray>

#include "pulsar.hh"

/**
 * Generate the reduced basis set by using a greedy algorithm.
 *
 * @param indices The indices of the pulsars in the data series
 * @param schedule The time values in the data series
 * @param pulsars The pulsar grid
 * @param params The parameter space. A 1D array of vectors in which there is a
 *    'linspace' type point array.
 * @param covarianceMatrix The covariance matrix to use for the inner product
 * @param epsilon The maximum error for the reduced basis
 * @param RB_out The reduced basis output array
 */
void greedyReducedBasis (std::vector<unsigned short> & indices,
                         std::vector<double> & schedule,
                         std::vector<Pulsar> & pulsars,
                         std::vector< std::vector<double> > & params,
                         std::vector<double> & covarianceMatrix,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_out);

/**
 * Translate an integer to a point in the parameter space. We do not store a grid, which
 * would  use lots of memory, we just return the parameter from the space via this
 * function
 *
 * @param idx This is an index in the parameter space (ranges from 0 to the total number
 *    of points in the parameter space
 * @param params The list of parameter values. It is basically a vector of ranges.
 * @param param_out The output array
 */
void idToParam (unsigned long idx, 
        std::vector<std::vector<double> > & params, 
        std::vector<std::vector<double> > & param_out);
