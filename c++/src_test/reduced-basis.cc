#include <vector>
#include <iostream>

#include "../src/reduced-basis.hh"
#include "../src/linalg.hh"

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
    std::vector<std::vector<double> > RB, RB_param;

    // A unit matrix and the error array
    std::vector<double> A = {1,0,0,0,1,0,0,0,1}, sigma;

    // Initialize some vectors in 3D
    TestG::dataSpace.clear();
    TestG::dataSpace.push_back({1,0,0});
    TestG::dataSpace.push_back({2,0,0});
    TestG::dataSpace.push_back({1,3,0});
    TestG::dataSpace.push_back({0,0,10});
    TestG::dataSpace.push_back({0,3,0});
    TestG::dataSpace.push_back({0,3,7});

    int N = TestG::dataSpace.size();
    linspace(TestG::paramSpace.at(0), 0, 1, N);

    greedyReducedBasis (N, getTestData, A, 0.000000001, RB_param, RB, sigma, false);

    int r = 0;

    if (not (RB.size() == 4 or RB.size() == 3)) {
        r++;
    }

    return r;
}

void getTestData (unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out) {
    params_out.resize(1);
    params_out[0] = idx;

    data_out.resize(3);
    data_out[0] = TestG::dataSpace.at(idx)[0];
    data_out[1] = TestG::dataSpace.at(idx)[1];
    data_out[2] = TestG::dataSpace.at(idx)[2];
}
