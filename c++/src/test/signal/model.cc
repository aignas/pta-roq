#include <vector>
#include <iostream>
#include <cmath>

#include "../../libs/linalg.hh"

#include "../../libs/pulsar.hh"
#include "../../libs/signal/vectors.hh"
#include "../../libs/signal/model.hh"

int test_pulsarTerm () {
    Pulsar P (1e+22, 0.3, 0.5, 1e-9, 0, 0, 0, 0);
    dvec e = { 1, 1e+19, 0.1, 0.1, 0.1 },
         i = { 3, 2, 1.1e-7 },
         vOmega, pUV, a, A_p;
    double t = 3.5e+08, L, delta, p;

    UnitVectors v (i.at(0), i.at(1));

    vOmega = v.Omega();
    pUV = P.getUnitVector();

    L = P.getDistance();

    delta = L * (1 + dotProduct(vOmega, pUV));

    amplitude(e, a);
    basis (t, i, pUV, A_p, delta);
    p = pulsarTerm (t, e, i, L, pUV);

    int r = 0;

    if (fabs( p - dotProduct(a,A_p) ) > 1e-320) {
        r++;
    }
    
    return r;
}
