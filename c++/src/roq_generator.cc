
#include <iterator>
#include <iostream>
#include <fstream>

#include <vector>

#include "pulsar.hh"
#include "linalg.hh"
#include "signal-model.hh"
#include "signal-sampling.hh"

#include "roq.hh"

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
    unsigned int pulsarNumber = 18;
    // The range initializer using C++11
    std::vector<double> range {50, 100, 0, _M_PI, -_M_PI, _M_PI};

    // Define some variables for bigger timescales
    const double week = 7 * 3600 * 24,
                 year = 52 * week,

                 t_final = 5 * year,
                 dt_min = 2 * week,
                 dt_max = 2 * week;

    std::vector<double> t_init (pulsarNumber, 0);

    // Variables needed for RB generation
    std::vector<double> C_inv,
                        Grammian,
                        interpolationMatrix,
                        sigma,
                        EIM_points,
                        templateNorms,
                        data,
                        data_tilda;
    std::vector<long> EIM_indices;
    std::vector<std::vector<double> > RB_params, RB;

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Prepare the pulsars" << std::endl;
    pulsarGrid::randomizeData (G::pulsars, pulsarNumber, range, 1e-15);
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, G::indices, G::Times);

    // Construct the parameter space
    G::params.resize(8);
    linspace(G::params[0], 1,1,1);
    linspace(G::params[1], 1e19, 1e19, 1);
    linspace(G::params[2], 0.1, 0.3, 1);
    linspace(G::params[3], 0.1, 0.3, 1);
    linspace(G::params[4], 0.1, 0.3, 1);
    linspace(G::params[5], 0.1, 3, 10);
    linspace(G::params[6], -3, 3, 10);
    linspace(G::params[7], 1.0e-9, 1.0e-8, 30);

    // Generate some mock data, so that we can construct the roq rule
    // Construct a parameter
    std::vector<double> source_one = {1,1e19,0.1,0.1,0.1,2.5,0.2,5e-9};
    std::vector<std::vector<double> > sources;
    sources.push_back(source_one);
    generateSample (data, G::pulsars, G::indices, G::Times, sources, true);

    // Construct the dimensionalities vector
    unsigned long totalNumber = 1;
    for (unsigned int i = 0; i < G::params.size(); i++) {
        G::dimensionalities.push_back(G::params.at(i).size());
        totalNumber *= G::dimensionalities.at(i);
    }

    // Set the error
    const double epsilon = 1e-35;

    try {
        std::cout << "Generate the covariance matrix" << std::endl;
        // Constructing a covariance matrix
        genCovarianceMatrix (C_inv, G::pulsars, G::indices, G::Times, true, false, false, false);

        std::cout << "Generate the RB" << std::endl;
        greedyReducedBasis (totalNumber, getData, C_inv, epsilon, RB_params, RB, Grammian, templateNorms, sigma, true);

        std::cout << "Generate the interpolation points: ";
        greedyEIMpoints (RB_params, RB, EIM_indices, EIM_points, interpolationMatrix);
        std::cout << "DONE" << std::endl;
        
        std::cout << "Generate the ROQ rule: ";
        constructROQ (data, data_tilda, EIM_indices, Grammian, interpolationMatrix);
        std::cout << "DONE" << std::endl;

        // Output the vectors to std::cout
        std::cout << "N,error,p1,p2,p3,p4,p5,p6,p7,p8" << std::endl;
        for (unsigned i = 0; i < sigma.size(); i++) {
            std::cout << i+1 << ", " << sigma.at(i);
            for (unsigned j = 0; j < RB_params[0].size(); j++) {
                std::cout << ", " << RB_params.at(i).at(j);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "N,eimid" << std::endl;
        for (unsigned i = 0; i < EIM_indices.size(); i++) {
            std::cout << i+1 << ", " << EIM_indices.at(i);
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "N, templateNorms" << std::endl;
        for (unsigned i = 0; i < templateNorms.size(); i++) {
            std::cout << i+1 << ", " << templateNorms.at(i);
            std::cout << std::endl;
        }
        std::cout << std::endl;
        
        std::cout << "signal, signal_tilda" << std::endl;
        for (unsigned i = 0; i < data_tilda.size(); i++) {
            std::cout << i+1 << ", " << data_tilda.at(i);
            std::cout << std::endl;
        }
        std::cout << std::endl;

    } catch (const char* msg) {
        // Add some colours
        std::cerr << msg << std::endl;
    }

    // Lets output the generated basis parameters and the error
    // FIXME write functions which output to a file

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
