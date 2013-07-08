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
double innerProduct (dvec& x, dvec& A, dvec& y);

/**
 * Construct a grammian
 *
 * @param G An output array
 * @param set A set of basis vectors
 * @param A A matrix to evaluate the inner product with
 */
void constructGrammian (dvec &G, std::vector<dvec>& set, dvec& A);

/**
 * Project a vector onto basis vectors
 *
 * @param projectee The vector to project
 * @param set The set of vector on which to project
 * @param A The matrix for the inner product
 * @param G_inv The inverse Grammian matrix for orthogonalization of the set
 * @params coeffs The coefficients for speeding up the calculation. If there are
 *   projection coefficients already calculated, then we can use them to speed up the
 *   calculation
 */
void projectionResidual (dvec & projectee, std::vector<dvec> & set, 
                         dvec & A, dvec & G_inv, 
                         dvec & coeffs);

/**
 * Norm of a vector
 *
 * @param x A vector to norm
 * @param A A matrix to norm it
 */
inline double norm (dvec& x, dvec& A) {
    return innerProduct(x, A, x);
}

/**
 * Matlab like linalg function
 *
 * @param array_out The output array (passed as a reference)
 * @param min A double representing a min value
 * @param max A double representing a max value
 * @param N Number of points in the array
 */
void linspace (dvec & array_out, double min, double max, const unsigned int N);

/**
 * Invert a matrix
 *
 * @param A a matrix with Row Major arrangement
 * @param The dimensionality of the matrix (it is assumed we have a square matrix)
 */
void inverse (dvec & A);
