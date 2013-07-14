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

void axpyProduct(double a, dvec & x, dvec & y) {
    cblas_daxpy(y.size(), a, &x[0], 1, &y[0], 1);
}

void matrixVectorProduct (dvec& A, dvec& x, dvec& y) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, x.size(), x.size(), 1, &A[0], x.size(), &x[0], 1, 0, &y[0], 1);
}

double innerProduct(dvec& x, dvec& A, dvec& y) {
    dvec tmp (x.size());

    matrixVectorProduct (A, y, tmp);

    double r = cblas_ddot(x.size(), &x[0], 1, &tmp[0], 1);

    return r;
}

void extendGrammianOptimized (dvec &G, dvec& G_inv, std::vector<dvec>& set, std::vector<dvec>& set_hat) {
    // Assign the old matrix to the new one
    unsigned int N = set.size();
    G_inv.resize(N*N);

    double tmp = 0;
    // We are extending the grammian in a clever way
    for (unsigned int i = 0; i < N - 1; i++) {
        tmp = dotProduct(set[N-1], set_hat[i]);
        G.emplace(G.begin() + i*N + (N-1), tmp);
        G.push_back(tmp);
    }

    tmp = dotProduct(set[N-1], set_hat[N-1]);
    G.push_back(tmp);

    G_inv = G;

    inverse(G_inv);
}

void constructGrammianOptimized (dvec &G, std::vector<dvec>& set, std::vector<dvec>& set_hat) {
    unsigned int N = set.size();
    G.resize(N*N);

    if (N != set_hat.size()) {
        std::cout << "The two sets should have the same number of members" << std::endl;
    }

    for (unsigned int i = 0; i < N; i++) {
        // The matrix is symmetric, hence calculate the diagonal terms here
        G.at(i * (N + 1)) = dotProduct (set[i], set_hat[i]);

        for (unsigned int j = i+1; j < N; j++) {
            // calculate the offdiagonal terms
            G.at(N*i + j) = dotProduct (set[i], set_hat[j]);
            G.at(i + N*j) = G[j + N*i];
        }
    }
}

void constructGrammian (dvec &G, std::vector<dvec>& set, dvec& A) {
    unsigned int N = set.size();
    G.resize(N*N);

    std::vector<dvec> set_tmp = set;

    for (unsigned int i = 0; i < N; i++) {
        matrixVectorProduct (A, set.at(i), set_tmp.at(i));
    }

    // Use the slightly optimized route
    constructGrammianOptimized (G, set, set_tmp);
}

void projectionResidual (dvec & projectee, std::vector<dvec> & set, 
                         dvec & A, dvec & G_inv) {
    unsigned int N = set.size();

    for (unsigned int i = 0; i < N; i++) {
        // Calculate the last coefficient, as we already are storring the rest
        double c_i = 0;
        for (unsigned int j = 0; j < N; j++) {
            c_i += G_inv[j + i*N] * innerProduct(projectee, A, set[j]);
        }

        axpyProduct(- c_i, set[i], projectee);
    }

}

void projectionResidualOptimized (dvec & projectee, std::vector<dvec> & set, 
                         std::vector<dvec> & set_hat, dvec & G_inv) {
    unsigned int N = set.size();

    double c_i;

    for (unsigned int i = 0; i < N; i++) {
        c_i = 0;
        for (unsigned int j = 0; j < N; j++) {
            // FIXME check for overflow/underflow
            c_i += G_inv[j + i*N] * dotProduct(projectee, set_hat[j]);
        }

        // FIXME check for overflow/underflow
        axpyProduct(- c_i, set[i], projectee);
    }

}

double projectionErrorStable (double projecteeNorm, dvec& projectee,
                            std::vector<dvec> & set_hat,
                            dvec & G, dvec & G_inv) {
    unsigned int N = set.size();

    dvec coef (N, 0);
    dvec projections (N, 0);

    // Parallelizable
    for (unsigned int i = 0; i < N; i++) {
        projections[i] = dotProduct(projectee, set_hat[i]);
    }

    // Calculate the coefficients, use BLAS
    matrixVectorProduct(G_inv, projection, coef);

    double error = projecteeNorm;

    error -= 2 * dotProduct(coef, projections);
    
    error += innerProduct (coef, G, coef);

    return error;
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
