#include <vector>
#include <cav/Vec3.hh>

class UnitVectors {
    private:
        double mTheta, mPhi;
    public:
        UnitVectors (double theta, double phi);
        UnitVectors~ ();

        cav::Vec3 Omega ();
        cav::Vec3 m ();
        cav::Vec3 n ();
}

UnitVectors::UnitVectors (double theta, double phi) {
    mTheta = theta;
    mPhi = phi;
}

UnitVectors::Omega () {
    return cav::Vec3 (
            -sin(mTheta)*cos(mPhi),
            -sin(mTheta)*sin(mPhi),
            -cos(mTheta))
}

UnitVectors::m () {
    return cav::Vec3 (
            -sin(mPhi),
            cos(mPhi),
            0)
}

UnitVectors::n () {
    return cav::Vec3 (
            -cos(mTheta)*cos(mPhi),
            -cos(mTheta)*sin(mPhi),
             sin(mTheta)
            )
}
