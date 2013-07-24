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

    if (argc != 7) {
        std::cerr << "Not enough parameters!!! The usage is as follows:\n"
                  << "\t " << argv[0] << " rc date N t_f dt_min dt_max\n\n"
                  << "Meanings of the options are:\n"
                  << "\t rc The configuration file\n"
                  << "\t date Time stamp in whatever format you want, but it should be preferably YYYY-MM-DD-HH-MM-SS\n"
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
    const unsigned int pulsarNumber = helper::convertToUnsignedInt(argv_s[2]);
    const double t_final = year * helper::convertToDouble(argv_s[3]),
                 dt_min = week * helper::convertToDouble(argv_s[4]),
                 dt_max = week * helper::convertToDouble(argv_s[5]);

    // Declare all the data structures
    // Randomize the starting times?
    std::vector<double> t_init (pulsarNumber, 0),
                        range = {1e20, 1e24, 0, 3.1, -3.1, 3.1},
                        Times, r;
    std::vector<unsigned short> indices;
    std::vector<Pulsar> pulsars;
    std::vector<std::vector<double> > sources;

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Generating a pulsar grid and the schedule for the measurements: "; std::cout.flush();
    pulsarGrid::randomizeData (pulsars, pulsarNumber, range, 5e-9);
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, indices, Times);
    std::cout << "DONE" << std::endl;

    // FIXME must be configurable from a file
    // Add a single source
    csv2arrayArrayDouble(fnames.at(0), sources, delim);

    // Generate a residual vector
    std::cout << "Generating a residual vector: "; std::cout.flush();
    try {
        generateSample (r, pulsars, indices, Times, sources);
        std::cout << "DONE" << std::endl;
    } catch (const std::bad_alloc&) {
        std::cerr << "Something went wrong" << std::endl;
        return 1;
    }

    // Output the data into files:
    // Output the pulsar parameters
    pulsar2csv (fnames.at(1), pulsars, delim);

    // Output the schedule
    arraysShortDouble2csv (fnames.at(2), indices, Times, delim);

    // Residuals out
    arrayDouble2csv (fnames.at(3), r, delim);

    return 0;
}
