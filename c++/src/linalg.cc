// This needs guards, as otherwise the compiler doesn't know what to do. More on it
// here:
// http://www.lindonslog.com/programming/atlas-blas-lapack-linear-algebra-libraries/
extern "C" {
#include <cblas.h>
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int*
            lwork, int* INFO);
}

// Local linear algebra routines to abstract things
#include "linalg.hh"
#include "signal-sampling.hh"

#include <iostream>

double dotProduct (dvec& x, dvec& y) {
    double r;
    if (x.size() == y.size()) {
        r = cblas_ddot(x.size(), &x[0], 1, &y[0], 1);
    } else {
        r = 0;
    }

    return r;
}

double innerProduct(dvec& x, dvec& A, dvec& y) {
    std::vector<double> tmp (x.size(), 0);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, y.size(), y.size(), 1, &A[0], y.size(), &y[0], 1, 0, &tmp[0], 1);

    double r = cblas_ddot(x.size(), &x[0], 1, &tmp[0], 1);

    return r;
}

void constructGrammian (dvec &G, std::vector<dvec>& set, dvec& A) {
    unsigned int N = set.size();
    G.resize(N*N);

    for (unsigned int i = 0; i < N; i++) {
        // The matrix is symmetric, hence calculate the diagonal terms here
        G.at(i * (N + 1)) = innerProduct (set[i], A, set[i]);

        for (unsigned int j = i+1; j < N; j++) {
            // calculate the offdiagonal terms
            G.at(N*i + j) = innerProduct (set[i], A, set[j]);
            G.at(i + N*j) = G[j + N*i];
        }
    }
}

void projectionResidual (dvec & projectee, std::vector<dvec> & set, 
                         dvec & A, dvec & G_inv, 
                         dvec & coeffs) {
    std::vector<double> tmp (projectee.size());
    unsigned int N = set.size(), Nc = coeffs.size();

    /*
    // Calculate the missing coefficients
    for (unsigned int i = Nc; i < N; i++) {
        coeffs.push_back(0);

        // Calculate the last coefficient, as we already are storring the rest
        for (unsigned int j = 0; j < N; j++) {
            coeffs[i] += G_inv[j + i*N] * innerProduct(projectee, A, set[j]);
        }
    }
    */

    for (unsigned int i = 0; i < N; i++) {
        // Calculate the last coefficient, as we already are storring the rest
        double c_i = 0;
        for (unsigned int j = 0; j < N; j++) {
            c_i += G_inv[j + i*N] * innerProduct(projectee, A, set[j]);
        }

        // Use BLAS to multiply by a constant
        cblas_daxpy(projectee.size(), - c_i, &set[i][0], 1, &projectee[0], 1);
    }

}

void linspace (std::vector<double> & array_out, double min, double max, const unsigned int N) {
    array_out.resize(N);

    // Check if the min and max are not mixed
    if (min > max) {
        double tmp = min;
        min = max;
        max = tmp;
    }

    // Push back the first element
    array_out[0] = min;

    if (N != 1) {
        const double dx = (max - min) / (N - 1);

        // push back the intermediate elements
        for (unsigned int i = 1; i < N - 1; i++) {
            array_out.at(i) = min + i*dx;
        }

        // Push back the last element
        array_out.at(N-1) = max;
    }
}

void inverse (std::vector<double> & A) {
    int N = sqrt(A.size());

    int INFO;

    std::vector<int> IPIV (N+1);
    int LWORK = A.size();
    std::vector<double> WORK (LWORK);

    // Execute the LAPACK routines
    dgetrf_(&N, &N, &A[0], &N, &IPIV[0], &INFO);
    dgetri_(&N, &A[0], &N, &IPIV[0], &WORK[0], &LWORK, &INFO);

    IPIV.clear();
    WORK.clear();
}

int arrayEqual (dvec &x, dvec &y) {
    int r = 0;
    bool equal = (x.size() == y.size());

    for (unsigned int i = 0; i < x.size() and equal; i++) {
        equal = equal and (x[i] - y[i] < 1e-300);
    }

    if (not equal) {
        r++;
    }

    return r;
}

void findMax (dvec & A, long & argmax_out, double & max_out) {
    // FIXME exception if the array is 0 size
    unsigned long N = A.size();

    // No need to compare first element with the first
    max_out = A.at(0), argmax_out = 0;

    /*
     * Parallelize the maximum search, taken from:
     * http://openmp.org/forum/viewtopic.php?f=3&t=70&start=0&hilit=max+reduction+array
     */
#pragma omp parallel for
    for (unsigned long i = 1; i < N; i++) {
        if (max_out < A.at(i)) {
#pragma omp critical
            {
                if (max_out < A.at(i)) {
                    argmax_out = i;
                    max_out = A.at(i);
                }
            }
        }
    }
}
