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
 * A matrix-vector product y = Ax
 *
 * @param A The matrix
 * @param x The multiplied vector
 * @param y The output vector
 */
void matrixVectorProduct (dvec& A, dvec& x, dvec& y);

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
 * An optimized version of constructing the grammian
 *
 * @param G An output array
 * @param set A set of basis vectors
 * @param set_hat A set of basis vectors premultiplied with the innerproduct matrix
 */
void constructGrammianOptimized (dvec &G, std::vector<dvec>& set, std::vector<dvec>& set_hat);

/**
 * Extend the Grammian matrix
 *
 * @param G The grammian to extend
 * @param set A set of basis vectors
 * @param set_hat A set of basis vectors premultiplied with the innerproduct matrix
 */
void extendGrammianOptimized (dvec &G, dvec &G_inv, std::vector<dvec>& set, std::vector<dvec>& set_hat);

/**
 * Project a vector onto basis vectors
 *
 * @param projectee The vector to project
 * @param set The set of vector on which to project
 * @param A The matrix for the inner product
 * @param G_inv The inverse Grammian matrix for orthogonalization of the set
 */
void projectionResidual (dvec & projectee, std::vector<dvec> & set, 
                         dvec & A, dvec & G_inv);

/**
 * Project a vector onto basis vectors, an optimized version, which uses to sets of
 * vectors instead of one set and a matrix
 *
 * @param projectee The vector to project
 * @param set The set of vector on which to project
 * @param set_hat The set of vector on which to project, multiplied by the innerproduct
 *   matrix.
 * @param G_inv The inverse Grammian matrix for orthogonalization of the set
 */
void projectionResidualOptimized (dvec & projectee, std::vector<dvec> & set, 
                         std::vector<dvec> & set_hat, dvec & G_inv);

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
 * @param A_inv The output (i.e. container for the inverse);
 */
void inverse (dvec & A);

/**
 * Check if arrays containing double are equal
 *
 * @param x Array
 * @param y Array
 */
int arrayEqual (dvec &x, dvec &y);

/**
 * Find a maximum value in the array
 *
 * @param A The array
 * @param argmax_out The index of the maximum value
 * @param max_out The value of the maximum
 */
void findMax (dvec & A, long & argmax_out, double & max_out);
