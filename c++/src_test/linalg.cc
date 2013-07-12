// Some general stuff
#include <vector>
#include <iostream>

// Local linear algebra routines to abstract things
#include "../src/linalg.hh"
#include "linalg.hh"

int test_dotProduct() {
    std::vector<double> x = {1,3,2}, y = {5,2,1};

    double p = dotProduct(x,y);
    int r = 0;

    // Check the answer
    if (p != 13) {
        r++;
    }

    return r;
}

int test_innerProduct () {
    std::vector<double> x = {1,2}, y = {-1, -3}, A = {1, 3, 2, 1};

    double p = innerProduct(x, A, y);
    int r = 0;

    // Check the answer
    if (p + 20 > 1e-50) {
        r++;
    }

    return r;
}

int test_inverse () {
    std::vector<double> A = {1, 2, 0, 4}, A_answer = {1, -0.5, 0, 0.25};

    inverse(A);

    return arrayEqual(A, A_answer);
}

int test_constructGrammian () {
    std::vector<double> x = {1,1,0}, y = {0,1,1}, z = {0,1,0}, 
                        A = {1,0,0,0,1,0,0,0,1}, G,
                        G_answer = {2,1,1,1,2,1,1,1,1};

    std::vector<std::vector<double> > set = {x,y,z};

    constructGrammian (G, set, A);

    return arrayEqual(G, G_answer);
}

int test_projectionResidual () {
    std::vector<double> x = {1,0.5,0}, y = {1,3,0}, z = {7,1,1},
                        A = {1,0,0,0,1,0,0,0,1}, G_inv, z_answer = {0,0,1};

    std::vector<std::vector<double> > set = {x,x,y};

    // Construct the inverse of the Grammian
    constructGrammian (G_inv, set, A);
    inverse(G_inv);

    // Check the projection
    projectionResidual (z, set, A, G_inv);

    return arrayEqual(z,z_answer);
}

int test_linspace () {
    std::vector<double> x, x_answer = {0,1,2,3,4};
    linspace(x, 0, 4, 5);

    return arrayEqual(x,x_answer);
}

int test_arrayEqual () {
    std::vector<double> x1 = {0,1,2}, x2 = {0,1,2}, y1 = {0,2,1};

    int r = 0;

    if (not (arrayEqual(x1,x2) == 0 and arrayEqual(x1,y1) == 1)) {
        r++;
    }

    return r;
}

int test_findMax () {
    // A semi random sequence
    std::vector<double> x = { 1,0,5,939,329348,8,98,98,97,5,3,3984,7398,938479,37894,384,739,749,9838794,7893978,4978,39874,9873897,3};

    long arg_max;
    double max;

    findMax(x,arg_max,max);

    int r = 0;

    if (not (arg_max == 22 and max - x[22] < 1e-300)) {
        r++;
    }

    return r;
}
