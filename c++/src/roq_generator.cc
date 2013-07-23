
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <locale>
#include <ctime>

#include <vector>

#include "libs/pulsar.hh"
#include "libs/linalg.hh"
#include "libs/signal/model.hh"
#include "libs/signal/sampling.hh"
#include "libs/iocsv.hh"
#include "libs/roq.hh"

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
    std::vector<std::string> fnames, argv_s;
    std::string delim;

    if (argc != 4) {
        std::cerr << "Not enough parameters!!! The usage is as follows:\n"
                  << "\t bin_name rc date date_in\n\n"
                  << "Meanings of the options are:\n"
                  << "\t rc The configuration file\n"
                  << "\t date Time stamp in whatever format you want, but it should be preferably YYYY-MM-DD-HH-MM-SS\n"
                  << "\t date_in Time stamp in whatever format you want. This should denote the time stamp for the file you want to read\n"
                  << std::endl;

        return 1;
    } else {
        for (int i = i; i < argc ; i++) {
            argv_s.push_back(argv[i]);
        }
    }

    parseRoqRC (argv_s[0], argv_s[1], argv_s[2], fnames, delim);

    // Variables needed for RB generation
    std::vector<double> C_inv, Grammian,
                        interpolationMatrix,
                        sigma,
                        EIM_points,
                        templateNorms,
                        data,
                        data_tilda,
                        Times_new;
    std::vector<long> EIM_indices;
    std::vector<unsigned short> indices_new;
    std::vector<std::vector<double> > RB_params, RB;

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Read pulsar and schedule data from a file: " << std::endl;
    csv2pulsar(fnames[0], G::pulsars, delim);

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

    // Read the schedule and the residuals from a file
    csv2arraysShortDouble(fnames[1], G::indices, G::Times, delim);
    csv2arrayDouble(fnames[2], data, delim);

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
        arrayArrayDouble2csv (fnames[3], RB_params, delim);
        arrayDouble2csv (fnames[4], sigma, delim);

        // Generate a new schedule with the eim points
        for (unsigned i = 0; i < EIM_indices.size(); i++) {
            indices_new.push_back(G::indices.at(EIM_indices.at(i)));
            Times_new.push_back(G::Times.at(EIM_indices.at(i)));
        }
        arraysShortDouble2csv (fnames[5], indices_new, Times_new, delim);
        arraysLongDouble2csv (fnames[6], EIM_indices, EIM_points, delim);

        arrayDouble2csv (fnames[7], data_tilda, delim);

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
