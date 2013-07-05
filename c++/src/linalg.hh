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
double dotProduct (dvec& x, dvec& y);

/**
 * A generalised inner product of 2 vectors of same dimensionality n, which are
 * multiplied with a given matrix.
 *
 * This function internally uses ATLAS, so should be fast enough
 *
 * @param x A vector of doubles
 * @param matrix a metrix of the same dimensionality as the vectors
 * @param y A vector of doubles
 */
double innerProduct (dvec& x, darray& A, dvec& y);

/**
 * Construct a grammian
 *
 * @param set A set of basis vectors
 * @param A A matrix to evaluate the inner product with
 */
darray constructGrammian (std::vector<dvec>& set, darray& A);

/**
 * Project a vector onto basis vectors
 *
 * @param projectee The vector to project
 * @param set The set of vector on which to project
 * @param A The matrix for the inner product
 * @param G The Grammian matrix for orthogonalization of the set
 */
dvec projectionIntoSet (dvec& projectee, std::vector<dvec>& set, darray& A, darray& G);

/**
 * Norm of a vector
 *
 * @param x A vector to norm
 * @param A A matrix to norm it
 */
inline double norm (dvec& x, darray& A);
