
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <locale>
#include <ctime>

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

/// 
// This will be a simple interface.
//
// The parameters, which can (should) be given are:
//      * timestamp: The timestamp of the data
//      * directory: The directory where to store the data
int main(int argc, char * argv []) {
    if (argc != 3) {
        std::cerr << "Not enough parameters. Please enter the timestamp and the saving directory for the data." << std::endl;

        return 1;
    }

    // Output filenames 
    std::vector<std::string> fnames = {
        "rb",       // Output reduced basis
        "time-roq"  // Output schedule (indices and Times)
        "data-roq", // Output roq rule
    };

    std::string prefix = "roq", 
                extension = ".csv";

    std::stringstream fnameFull;
    for (unsigned i = 0; i < fnames.size(); i++) {
        fnameFull.str(std::string());
        fnameFull << argv[2] << prefix << "-" << argv[1] << "-"  << fnames[i] << extension;
        fnames[i] = fnameFull.str();
    }
    std::cout << std::endl;

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

    // FIXME Read the pulsar data from a file
    // Randomize the pulsar structure and generate a schedule
    std::cout << "Read pulsar and schedule data from a file: " << std::endl;

    // FIXME check if it is possible to have reduced basis for the Basis functions of
    // the signal (A^i):
    G::params.resize(8);
    linspace(G::params[0], 1,1,1);
    linspace(G::params[1], 1e19, 1e19, 1);
    linspace(G::params[2], 0.1, 0.3, 1);
    linspace(G::params[3], 0.1, 0.3, 1);
    linspace(G::params[4], 0.1, 0.3, 1);
    linspace(G::params[5], 0.1, 3, 10);
    linspace(G::params[6], -3, 3, 10);
    linspace(G::params[7], 2e-9, 4e-7, 30);

    // Construct the dimensionalities vector
    unsigned long totalNumber = 1;
    for (unsigned int i = 0; i < G::params.size(); i++) {
        G::dimensionalities.push_back(G::params.at(i).size());
        totalNumber *= G::dimensionalities.at(i);
    }

    // FIXME Read some data from a file

    // Set the error
    const double epsilon = 1e-35;

    // Protected code
    try {
        std::ofstream fout;
        std::cout << "Generate the covariance matrix: "; std::cout.flush();
        // Constructing a covariance matrix
        genCovarianceMatrix (C_inv, G::pulsars, G::indices, G::Times, true, false, false, false);
        std::cout << "DONE" << std::endl;

        std::cout << "Generate the RB: "; std::cout.flush();
        greedyReducedBasis (totalNumber, getData, C_inv, epsilon, RB_params, RB, Grammian, templateNorms, sigma, true);
        std::cout << "DONE" << std::endl;

        std::cout << "Generate the interpolation points: "; std::cout.flush();
        greedyEIMpoints (RB_params, RB, EIM_indices, EIM_points, interpolationMatrix);
        std::cout << "DONE" << std::endl;
        
        std::cout << "Generate the ROQ rule: "; std::cout.flush();
        constructROQ (data, data_tilda, EIM_indices, Grammian, interpolationMatrix);
        std::cout << "DONE" << std::endl;

        // Output the basis params and error to a file
        fout.open(fnames[0]);
        fout << "error";
        for (unsigned j = 0; j < RB_params[0].size(); j++) {
            fout << ", p" << j;
        }
        fout << std::endl;
        for (unsigned i = 0; i < sigma.size(); i++) {
            fout << sigma.at(i);
            for (unsigned j = 0; j < RB_params[0].size(); j++) {
                fout << ", " << RB_params.at(i).at(j);
            }
            fout << std::endl;
        }
        fout.close();

        // Generate a new schedule with the eim points
        fout.open(fnames[1]);
        fout << "eimid, schedid, schedtime" << std::endl;
        for (unsigned i = 0; i < EIM_indices.size(); i++) {
            fout << EIM_indices.at(i) << ", " <<  
                    G::indices.at(EIM_indices.at(i)) << ", " << 
                    G::Times.at(EIM_indices.at(i)) << 
                    std::endl;
        }
        fout.close();

        fout.open(fnames[3]);
        fout << "signal_tilda" << std::endl;
        for (unsigned i = 0; i < data_tilda.size(); i++) {
            fout << data_tilda.at(i) << std::endl;
        }
        fout.close();

    } catch (const char* msg) {
        // Add some colours
        std::cerr << msg << std::endl;
    }

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
