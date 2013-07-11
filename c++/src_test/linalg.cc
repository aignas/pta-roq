// Some general stuff
#include <vector>
#include <valarray>
#include <iostream>

// Local linear algebra routines to abstract things
#include "../src/linalg.hh"

int test_dotProduct() {
    std::vector<double> x (3), y(3);
    x[0] = 1, y[0] = 5;
    x[1] = 3, y[1] = 2;
    x[2] = 2, y[2] = 1;

    double p = dotProduct(x,y);
    int r = 0;

    // Check the answer
    if (p != 13) {
        r++;
    }

    return r;
}

int test_innerProduct () {
    std::vector<double> x (2), y(2);
    x[0] = 1, y[0] = -1;
    x[1] = 2, y[1] = -3;
    std::valarray<double> A (4);
    A[0] = 1, A[1] = 3, A[2] = 2, A[3] = 1;

    double p = innerProduct(x, A, y);
    int r = 0;

    // Check the answer
    if (p != 20) {
        r++;
    }

    return r;
}


