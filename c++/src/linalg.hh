/**
 * This will be the linear algebra wrapper for the RB algorithms
 */

#include <vector>
#include <valarray>

typedef std::vector<double> dvec;
typedef std::valarray<double> darray;

/**
 * A generalised dot product of 2 vectors of same dimensionality n.
 *
 * This function internally uses ATLAS, so should be fast enough
 *
 * @param x A vector of doubles
 * @param y A vector of doubles
 */
double dotProduct (dvec x, dvec y);

/**
 * A generalised inner product of 2 vectors of same dimensionality n, which are
 * multiplied with a given matrix.
 *
 * This function internally uses ATLAS, so should be fast enough
 *
 * @param x A vector of doubles
 * @param y A vector of doubles
 * @param matrix a metrix of the same dimensionality as the vectors
 */
double innerProduct (dvec x, dvec y, darray matrix);
