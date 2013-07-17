#include <vector>
#include <iostream>

#include "roq.hh"

#include "../src/roq.hh"
#include "../src/linalg.hh"
#include "../src/random-helper.hh"

//////////////////////
// This is needed for the test_greedyReducedBasis
namespace TestG {
    std::vector<std::vector<double> > paramSpace (1);

    std::vector<std::vector<double> > dataSpace;
}

void getTestData (unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out);
//////////////////////

int test_idToList () {
    unsigned long N = 3*5*8;

    std::vector<unsigned int> dimensionalities = {3,5,8};
    std::vector<std::vector<unsigned int> > ids;
    std::vector<unsigned int> list (3);


    // Generate ids by hand
    for (unsigned int k = 0; k < 8; k++) {
        for (unsigned int j = 0; j < 5; j++) {
            for (unsigned int i = 0; i < 3; i++) {
                list = {i,j,k};
                ids.push_back(list);
            }
        }
    }

    bool match = true;

    for (unsigned int i = 0; i < N and match; i++) {
        list = {0, 0, 0};
        idToList (i, dimensionalities, list);
        for (unsigned int j = 0; j < 3; j++) {
            match = match and (ids[i][j] == list[j]);
        }
    }

    int r = 0;

    if (not match) {
        r++;
    }

    return r;
}

int test_greedyReducedBasis () {
    // The dimensions of the vector space
    unsigned N = 116;

    // Initialize some vectors in N-D
    TestG::dataSpace.clear();
    for (unsigned i = 0; i < N*5; i++) {
        std::vector<double> tmp (N);
        for (unsigned j = 0; j < N; j++) {
            tmp.at(j) = random_uniform (-1, 1);
        }
        TestG::dataSpace.push_back(tmp);
    }

    linspace(TestG::paramSpace.at(0), 0, 1, N*15);

    // Set the error
    const double epsilon = 1e-24;

    // Generate some symmetric matrix:
    std::vector<double> A (N*N);

    for (unsigned i = 0; i < N; i++) {
        //A.at(i*(N + 1)) = 1./random_uniform (1, 100);
        A.at(i*(N + 1)) = 1;
    }

    // Generate a residual vector
    std::vector<std::vector<double> > RB_params, RB;
    std::vector<double> sigma, G, tmp;
    greedyReducedBasis (TestG::dataSpace.size(), getTestData, A, epsilon, RB_params, RB, G, tmp, sigma, false);

    random_uniform_free();

    int r = 0;

    if (not (RB.size() == N or RB.size() == N + 1)) {
        r++;
    }

    return r;
}

void getTestData (unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out) {
    params_out.resize(1);
    params_out[0] = idx;

    data_out = TestG::dataSpace.at(idx);
}
