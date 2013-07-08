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
            G.at(N*j + i) = innerProduct (set[i], A, set[j]);
            G.at(j + N*i) = G[i + N*j];
        }
    }
}

void projectionResidual (dvec & projectee, std::vector<dvec> & set, 
                         dvec & A, dvec & G, 
                         dvec & coeffs) {
    std::vector<double> tmp (projectee.size());
    unsigned int N = set.size();

    // Calculate the missing coefficients
    for (unsigned int i = coeffs.size(); i < N; i++) {
        coeffs.push_back(0);

        // Calculate the last coefficient, as we already are storring the rest
        for (unsigned int j = 0; j < N - 1; j++) {
            coeffs.at(i) += G[i + j*N] * innerProduct(projectee, A, set[j]);
        }
    }

    for (unsigned int i = 0; i < N; i++) {
        // Use BLAS to multiply by a constant
        cblas_daxpy(projectee.size(), - coeffs[i], &set[i][0], 1, &projectee[0], 1);
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
