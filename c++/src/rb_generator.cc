
#include <vector>
#include <valarray>
#include <iostream>
#include <fstream>

#include "pulsar.hh"
#include "linalg.hh"
#include "signal-sampling.hh"

#include "reduced-basis.hh"

/*
 * Use global variables here and although most of the times the usage of them is not
 * justified, this way the code readability increases a lot and I also prefix them with
 * a namespace name. This should make it easy to understand weather a variable is global
 * or not
 */
namespace G {
    // The pulsar grid
    std::vector<Pulsar> pulsars;

    // Time schedule and indices of the pulsars for each data point
    std::vector<double> Times;
    std::vector<unsigned short> indices;

    // The dimensionalities of the parameter space
    std::vector<unsigned int> dimensionalities;

    // The parameterspace
    std::vector<std::vector<double> > params;
}

// Declare the get data function, which will be used in the reduced basis generator
void getData(unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out);

int main() {
    unsigned int pulsarNumber = 36;
    // The range initializer using C++11
    std::vector<double> range {50, 100, 0, _M_PI, -_M_PI, _M_PI};

    // Define some variables for bigger timescales
    const double week = 7 * 3600 * 24,
                 year = 52 * week,

                 t_final = 5 * year,
                 dt_min = 2 * week,
                 dt_max = 2 * week;

    std::vector<double> t_init (pulsarNumber, 0);

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Randomizing the initial conditions" << std::endl;
    pulsarGrid::randomizeData (G::pulsars, pulsarNumber, range, 1e-15);
    std::cout << "Constructing the schedule of measuremnts" << std::endl;
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, G::indices, G::Times);

    // Construct the parameter space
    G::params.resize(8);
    linspace(G::params[0], 1,1,1);
    linspace(G::params[1], 1e9, 1e9, 1);
    linspace(G::params[2], 0.1, 0.3, 1);
    linspace(G::params[3], 0.1, 0.3, 1);
    linspace(G::params[4], 0.1, 0.3, 1);
    linspace(G::params[5], 1.2, 1.5, 10);
    linspace(G::params[6], -0.2, 0.2, 10);
    linspace(G::params[7], 5.0e-9, 1.0e-8, 30);

    // Construct the dimensionalities vector
    unsigned long totalNumber = 1;
    for (unsigned int i = 0; i < G::params.size(); i++) {
        G::dimensionalities.push_back(G::params.at(i).size());
        totalNumber *= G::dimensionalities.at(i);
    }

    // Set the error
    const double epsilon = 1e-35;

    // Constructing a covariance matrix
    std::cout << "Constructing the covariance matrix" << std::endl;
    std::vector<double> C_inv;
    genCovarianceMatrix (C_inv, G::pulsars, G::indices, G::Times, true, false, false, false);

    // Generate a residual vector
    std::cout << "Generating the reduced basis set" << std::endl;
    std::vector<std::vector<double> > RB_params, RB;
    std::vector<double> sigma;
    greedyReducedBasis (totalNumber, getData, C_inv, epsilon, RB_params, RB, sigma, true);

    // Lets output the generated basis parameters and the error
    std::cout << "Outputing the generated basis to a files" << std::endl;
    std::ofstream sigma_out, RB_params_out;
    sigma_out.open("sigma.out");
    RB_params_out.open("params.out");
    for (unsigned int i = 0; i < sigma.size(); i++) {
        sigma_out << sigma.at(i) << std::endl;
        for(unsigned int j = 0; j < RB_params[i].size(); j++) {
            RB_params_out << RB_params[i][j] << " ";
        }
        RB_params_out << std::endl;
    }
    sigma_out.close();
    RB_params_out.close();

    return 0;
}

/**
 * This function will return data to the reduced basis algorithm
 */
void getData(unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out) {
    std::vector<unsigned int> list;
    idToList(idx, G::dimensionalities, list);
    params_out.resize(list.size());

    // Construct a parameter
    for (unsigned int i = 0; i < list.size(); i++) {
        params_out.at(i) = G::params.at(i).at(list.at(i));
    }

    std::vector<std::vector<double> > sources (1, params_out);

    // Do not include noise as we are generating a template
    generateSample (data_out, G::pulsars, G::indices, G::Times, sources, false);
}
