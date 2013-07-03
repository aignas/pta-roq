// This needs guards, as otherwise the compiler doesn't know what to do. More on it
// here:
// http://www.lindonslog.com/programming/atlas-blas-lapack-linear-algebra-libraries/
//extern "C" {
//#include <cblas.h>
//}

#include <cmath>
#include <signal-vectors.hh>

UnitVectors::UnitVectors (double theta, double phi) {
    mTheta = theta;
    mPhi = phi;
}

std::vector<double> UnitVectors::Omega () {
    std::vector<double>  x(3);
    x[0] = -sin(mTheta) * cos(mPhi);
    x[1] = -sin(mTheta) * sin(mPhi);
    x[2] = -cos(mTheta);
    return x;
}

std::vector<double> UnitVectors::m () {
    std::vector<double>  x(3);
    x[0] = -sin(mPhi);
    x[1] =  cos(mPhi);
    x[2] = 0;
    return x;
}

std::vector<double> UnitVectors::n () {
    std::vector<double>  x(3);
    x[0] = -cos(mTheta) * cos(mPhi);
    x[1] = -cos(mTheta) * sin(mPhi);
    x[2] = -sin(mTheta);
    return x;
}
