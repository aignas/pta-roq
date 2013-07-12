
#include <vector>
#include <valarray>
#include <iostream>
#include <fstream>

#include "../src/linalg.hh"
#include "../src/random-helper.hh"

#include "../src/reduced-basis.hh"

/*
 * Use global variables here and although most of the times the usage of them is not
 * justified, this way the code readability increases a lot and I also prefix them with
 * a namespace name. This should make it easy to understand weather a variable is global
 * or not
 */
namespace G {
    std::vector<std::vector<double> > paramSpace (1);

    std::vector<std::vector<double> > dataSpace;
}

// Declare the get data function, which will be used in the reduced basis generator
void getData(unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out);

int main() {
    // The dimensions of the vector space
    unsigned N = 15;

    // Initialize some vectors in N-D
    G::dataSpace.clear();
    for (unsigned i = 0; i < N*15; i++) {
        std::vector<double> tmp (N);
        for (unsigned j = 0; j < N; j++) {
            tmp.at(j) = random_uniform (-30, 30);
            std::cout << tmp.at(j) << " ";
        }
        G::dataSpace.push_back(tmp);
        std::cout << std::endl;
    }

    for (unsigned i = 0; i < 15; i++) {
        std::cout << random_uniform() << " ";
    }
    std::cout << std::endl;

    linspace(G::paramSpace.at(0), 0, 1, N*15);

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
    std::vector<double> sigma;
    greedyReducedBasis (G::dataSpace.size(), getData, A, epsilon, RB_params, RB, sigma, true);

    random_uniform_free();

    return 0;
}

void getData (unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out) {
    params_out.resize(1);
    params_out[0] = idx;

    data_out = G::dataSpace.at(idx);
}
