
#include <vector>
#include <valarray>

#include "pulsar.hh"

/**
 * Generate the reduced basis set by using a greedy algorithm. This does not assume on
 * any functions whilst calculating the 
 *
 * @param N The number of points in the parameter space
 * @param getData A function of int, param_out and data_out variable, which would abstract the way
 *  I am getting data fram the parameter space.
 * @param A The matrix to use in the inner product
 * @param epsilon The maximum error for the reduced basis
 * @param RB_out The reduced basis output array
 * @param sigma_out The array for the error variation during the recursion
 */
void greedyReducedBasis (const unsigned long N,
                         void (*getData)(unsigned long idx, 
                                         std::vector<double> & params_out, 
                                         std::vector<double> & data_out),
                         std::vector<double> & A,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_param_out,
                         std::vector<std::vector<double> > & RB_out,
                         std::vector<double> & sigma_out);

/**
 * Generate a list of numbers from a single number. This acts as a good way to recover
 * data from a parameter space stored in an array sequentially
 *
 * @param idx A long number of choice, which should be in a range of the parameter
 *   space.
 * @param dim A list of dimensionalities of the parameter space.
 * @param list_out The output array
 */
void idToList (unsigned long idx, 
               std::vector<unsigned int> dim, 
               std::vector<unsigned int> & list_out);
