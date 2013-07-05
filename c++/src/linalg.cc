// This needs guards, as otherwise the compiler doesn't know what to do. More on it
// here:
// http://www.lindonslog.com/programming/atlas-blas-lapack-linear-algebra-libraries/
extern "C" {
#include <cblas.h>
#include <clapack.h>
}

// Local linear algebra routines to abstract things
#include "linalg.hh"

double dotProduct (dvec& x, dvec& y) {
    double r;
    if (x.size() == y.size()) {
        r = cblas_ddot(x.size(), &x[0], 1, &y[0], 1);
    } else {
        r = 0;
    }

    return r;
}

double innerProduct(dvec& x, darray& A, dvec& y) {
    double r;

    // Fix it
    if (x.size() == y.size() and x.size()*y.size() == A.size()) {
        std::vector<double> tmp (x.size(), 0);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, x.size(), x.size(), 1, &A[0], 1, &y[0], 1, 1, &tmp[0], 1);
        r = cblas_ddot(x.size(), &x[0], 1, &tmp[0], 1);
    } else {
        r = 0;
    }

    return r;
}

darray constructGrammian (std::vector<dvec>& set, darray& A) {
    unsigned int N = set.size();

    std::valarray<double> G (N*N);

    for (unsigned int i = 0; i < N; i++) {
        // The matrix is symmetric, hence calculate the diagonal terms here
        G[i * (N + 1)] = innerProduct (set[i], A, set[i]);

        for (unsigned int j = i+1; j < N; j++) {
            // calculate the offdiagonal terms
            G[N*j + i] = innerProduct (set[i], A, set[j]);
            G[j + N*i] = G[i + N*j];
        }
    }

    return G;
}

dvec projectionIntoSet (dvec& projectee, std::vector<dvec>& set, darray& A, darray& G) {
    unsigned int N = set.size();
    std::vector<double> tmp (projectee.size(), 0);
    double c_i;

    for (unsigned int i = 0; i < N; i++) {
        c_i = 0;
        for (unsigned int j = 0; j < N; j++) {
            c_i += G[i + j*N] * innerProduct(projectee, A, set[j]);
        }
        cblas_daxpy(projectee.size(), c_i, &set[i][0], 1, &tmp[0], 1);
    }

    return tmp;
}

inline double norm (dvec& x, darray& A) {
    return innerProduct(x, A, x);
}
