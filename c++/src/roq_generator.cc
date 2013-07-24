
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
    std::vector<std::string> fnames, argvs;
    std::string delim;
    for (int i = 0; i < argc ; i++) {
        argvs.push_back(argv[i]);
    }

    if (argc != 6) {
        std::cerr << "Not enough parameters!!! The usage is as follows:\n"
                  << "\t " << argvs[0] << " rc param_range_rc date date_in error\n\n"
                  << "Meanings of the options are:\n"
                  << "\t rc The configuration file\n"
                  << "\t param_range_rc The configuration file for the parameter ranges\n"
                  << "\t date Time stamp in whatever format you want, but it should be preferably YYYY-MM-DD-HH-MM-SS\n"
                  << "\t date_in Time stamp in whatever format you want. This should denote the time stamp for the file you want to read\n"
                  << "\t error An error to achieve in the RB approximation\n"
                  << std::endl;

        return 1;
    } else {
        std::cout << "\nStarting the ROQ construction" << std::endl;
    }

    if (parseRoqRC (argvs[1], argvs[3], argvs[4], fnames, delim) == 1) {
        std::cerr <<        "The syntax of the config file was incorrect. The syntax is described bellow: "
            << std::endl << "The key (i.e. the first value) matters here! But the whitespace doesn't unless it is a"
            << std::endl << "commented line. Also, the order in which the options are given doesn't matter either."
            << std::endl << "Commented lines start either with #."
            << std::endl << ""
            << std::endl << "Available keys:"
            << std::endl << "\t- ext The extension of the files"
            << std::endl << "\t- delim The delimiter used"
            << std::endl << "\t- sep The separator used in the filename"
            << std::endl << ""
            << std::endl << "\t- dir.in The output directory"
            << std::endl << "\t- prefix.in The prefix for the data input files"
            << std::endl << "\t- infile.pulsar Input pulsar parameters"
            << std::endl << "\t- infile.sched Input time stamps and indices of the timed pulsars"
            << std::endl << "\t- infile.resid Input the residuals"
            << std::endl << ""
            << std::endl << "\t- dir.out The output directory"
            << std::endl << "\t- prefix.out The prefix for the data output files"
            << std::endl << "\t- outfile.sched Ouput the modified schedule"
            << std::endl << "\t- outfile.eim Output the EIM points and indices"
            << std::endl << "\t- outfile.resid Ouput the quadrature rule"
            << std::endl << "\t- outfile.rb Ouput the reduced basis"
            << std::endl << "\t- outfile.rbparams Ouput the parameters used to construct the reduced basis"
            << std::endl;
    }

    // Set the error
    const double epsilon = helper::convertToDouble(argvs[5]);

    // Variables needed for RB generation
    std::vector<double> C_inv, Grammian, interpolationMatrix,
                        sigma, EIM_points, templateNorms,
                        data, data_tilda,
                        Times_new,
                        params_min, params_max;
    std::vector<long> EIM_indices;
    std::vector<unsigned short> indices_new;
    std::vector<unsigned int> params_N;
    std::vector<std::vector<double> > RB_params, RB;

    // Read pulsar and schedule data from a file
    csv2arraysShortDouble(fnames[1], G::indices, G::Times, delim);
    csv2pulsar(fnames[0], G::pulsars, delim);

    // FIXME check if it is possible to have reduced basis for the Basis functions of
    // the signal (A^i):
    csv2paramRanges(argvs[2], params_min, params_max, params_N, delim);
    G::params.resize(params_min.size());
    std::cout << "Generating the parameter space: " << std::endl;
    for (unsigned i = 0; i < G::params.size(); i++) {
        std::cout << "\t" << i << " " << params_N[i] << " points in range [" << params_min[i] << "; " << params_max[i] << "]." << std::endl;
        linspace(G::params[i], params_min[i], params_max[i], params_N[i]);
    }

    // Construct the dimensionalities vector
    unsigned long totalNumber = 1;
    for (unsigned int i = 0; i < G::params.size(); i++) {
        G::dimensionalities.push_back(G::params.at(i).size());
        totalNumber *= G::dimensionalities.at(i);
    }

    // Read the residuals from a file
    csv2arrayDouble(fnames[2], data, delim);

    // Protected code
    try {
        genCovarianceMatrix (C_inv, G::pulsars, G::indices, G::Times, true, false, false, false);
        greedyReducedBasis (totalNumber, getData, C_inv, epsilon, RB_params, RB, Grammian, templateNorms, sigma, true);
        greedyEIMpoints (RB_params, RB, EIM_indices, EIM_points, interpolationMatrix);
        constructROQ (data, data_tilda, EIM_indices, Grammian, interpolationMatrix);

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
