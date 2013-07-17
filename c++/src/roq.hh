#include <vector>

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
 * @param verbose should data be printed to stdout?
 */
void greedyReducedBasis (const unsigned long N,
                         void (*getData)(unsigned long idx, 
                                         std::vector<double> & params_out, 
                                         std::vector<double> & data_out),
                         std::vector<double> & A,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_param_out,
                         std::vector<std::vector<double> > & RB_out,
                         std::vector<double> & sigma_out,
                         bool verbose = false);

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
               std::vector<unsigned int> & dim, 
               std::vector<unsigned int> & list_out);

/**
 * Generate the reduced basis set by using a greedy algorithm. This does not assume on
 * any functions whilst calculating the 
 *
 * @param RB_param The reduced basis parameter array
 * @param RB The reduced basis array
 * @param indices_out The output array for the EIM indices
 * @param points_out The output array for the EIM points
 */
void greedyInterpolant (std::vector<std::vector<double> > & RB_param,
                        std::vector<std::vector<double> > & RB,
                        std::vector<long> & indices_out,
                        std::vector<double> points_out);

/**
 * Generate a quadrature rule vector
 *
 * @param r The input vector, which will be overwritten with the quadrature rule vector
 * @param RB The reduced basis array
 * @param idx_EIM The Empirical interpolation points
 */
void constructROQ (std::vector<double> r,
                   std::vector<std::vector<double> > & RB,
                   std::vector<long> & indices_out);
