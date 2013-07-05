
#include <vector>
#include <valarray>

// This needs guards, as otherwise the compiler doesn't know what to do. More on it
// here:
// http://www.lindonslog.com/programming/atlas-blas-lapack-linear-algebra-libraries/
extern "C" {
#include <cblas.h>
#include <clapack.h>
}

// Local linear algebra routines to abstract things
#include "linalg.hh"

typedef std::vector<double> dvec;
typedef std::valarray<double> darray;

double dotProduct (dvec x, dvec y) {
    double r;
    if (x.size() == y.size()) {
        r = cblas_ddot(x.size(), &x[0], 1, &y[0], 1);
    } else {
        r = 0;
    }

    return r;
}

double innerProduct(dvec x, dvec y, darray A) {
    double r;

    // Fix it
    if (x.size() == y.size()) {
        r = cblas_ddot(x.size(), &x[0], 1, &y[0], 1);
    } else {
        r = 0;
    }

    return r;
}
