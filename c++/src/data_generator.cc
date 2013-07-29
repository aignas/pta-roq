#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "libs/signal/model.hh"
#include "libs/signal/sampling.hh"
#include "libs/pulsar.hh"
#include "libs/iocsv.hh"

int main(int argc, char * argv []) {
    std::vector<std::string> fnames, argv_s;
    std::string delim;

    // Define some variables for bigger timescales
    // These should be in python or something similar.
    const double week = 7 * 3600 * 24,
                 year = 52 * week;

    if (argc != 8) {
        std::cerr << "Not enough parameters!!! The usage is as follows:\n"
                  << "\t " << argv[0] << " rc date rc_source N t_f dt_min dt_max\n\n"
                  << "Meanings of the options are:\n"
                  << "\t rc The configuration file\n"
                  << "\t date Time stamp in whatever format you want, but it should be preferably YYYY-MM-DD-HH-MM-SS\n"
                  << "\t rc_source The parameters for the source\n"
                  << "\t N The pulsar number for the simulation\n"
                  << "\t t_f Final simulation time in years\n"
                  << "\t dt_min Minimum interval between measurements for some pulsar (in weeks)\n"
                  << "\t dt_max Maximum interval between measurements for some pulsar (in weeks)\n"
                  << "\t The time schedule can be a bit randomized if dt_min != dt_max"
                  << std::endl;

        return 1;
    } else {
        for (int i = 1; i < argc ; i++) {
            argv_s.push_back(argv[i]);
        }
    }

    parseDataRC (argv_s[0], argv_s[1], fnames, delim);

    // Read the parameters
    const unsigned int pulsarNumber = helper::convertToUnsignedInt(argv_s[3]);
    const double t_final = year * helper::convertToDouble(argv_s[4]),
                 dt_min = week * helper::convertToDouble(argv_s[5]),
                 dt_max = week * helper::convertToDouble(argv_s[6]);

    // Declare all the data structures
    // Randomize the starting times?
    std::vector<double> t_init (pulsarNumber, 0),
    // FIXME must be configurable from a file
                        range = {1e+10, 1e+16, 0.01, 3.1, -3.1, 3.1},
                        Times, r;
    std::vector<unsigned short> indices;
    std::vector<Pulsar> pulsars;
    std::vector<std::vector<double> > sources;

    // Randomize the pulsar structure and generate a schedule
    pulsarGrid::randomizeData (pulsars, pulsarNumber, range, 5e-9);
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, indices, Times);

    // Add a single source
    csv2arrayArrayDouble(argv_s[2], sources, delim);

    // Generate a residual vector
    try {
        for (unsigned i = 0; i < sources.size(); i++) {
            if (i==0) {
                std::cout << "Generating the residuals with parameters:";
            } else {
                std::cout << "                                         ";
            }
            for (unsigned j = 0; j < sources[i].size(); j++) {
                std::cout << " " << sources[i][j];
            }
            std::cout << std::endl;
        }
        generateSample (r, pulsars, indices, Times, sources);
    } catch (const std::bad_alloc&) {
        std::cerr << "Something went wrong" << std::endl;
        return 1;
    }

    // Output the data into files:
    // Output the pulsar parameters
    pulsar2csv (fnames.at(0), pulsars, delim);
    // Output the schedule
    arraysShortDouble2csv (fnames.at(1), indices, Times, delim);
    // Residuals out
    arrayDouble2csv (fnames.at(2), r, delim);

    return 0;
}
