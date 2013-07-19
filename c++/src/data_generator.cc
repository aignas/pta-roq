#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "signal-model.hh"
#include "signal-sampling.hh"
#include "pulsar.hh"

int main(int argc, char * argv []) {
    // argv[1] = saving dir
    // argv[2] = timestamp
    if (argc != 3) {
        std::cerr << "Not enough parameters. Please enter the timestamp and the saving directory for the data." << std::endl;

        return 1;
    }

    const unsigned int pulsarNumber = 36;

    // Output filenames 
    std::vector<std::string> fnames = {
        "params",   // Output schedule (indices and Times)
        "pulsars",  // Output pulsar properties
        "schedule"  // Output schedule (indices and Times)
    };

    std::string prefix = "model", 
                separator = "-",
                extension = ".csv";

    std::stringstream fnameFull;
    // Generate the filenames
    for (unsigned i = 0; i < fnames.size(); i++) {
        fnameFull.str(std::string());
        fnameFull << argv[1] << "/" << prefix << separator << argv[2] << separator << fnames[i] << extension;
        fnames[i] = fnameFull.str();
    }

    // Define some variables for bigger timescales
    const double week = 7 * 3600 * 24,
                 year = 52 * week,

                 t_final = 5 * year,
                 dt_min = 2 * week,
                 dt_max = 2 * week;

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

    // Add a single source
    sources.push_back({1,1e19, 0.1, 0.1, 0.1, 0.3, 0.5, 1e-7});

    // Generate a residual vector
    std::cout << "Generating a residual vector: "; std::cout.flush();
    try {
        generateSample (r, pulsars, indices, Times, sources);
        std::cout << "DONE" << std::endl;
    } catch (const std::bad_alloc&) {
        std::cerr << "Something went wrong" << std::endl;
        return 1;
    }

    // Output the data and the pulsar schedule into a file:
    //      * params - parameters used to generate the data
    //      * pulsars - Pulsar props
    //      * schedule-i - Schedule for the ith pulsar
    std::ofstream fout;
    fout.open(fnames[0]);
    for (unsigned i = 0; i < sources.size(); i++) {
        for (unsigned j = 0; j < sources.at(i).size(); j++) {
            fout << sources.at(i).at(j) << " ";
        }
        fout << std::endl;
    }
    fout.close();

    fout.open(fnames[1]);
    std::vector<double> tmp;
    for (unsigned i = 0; i < pulsars.size(); i++) {
        pulsars[i].getAll(tmp);
        for (unsigned j = 0; j < tmp.size(); j++) {
            fout << tmp.at(j) << " ";
        }
        fout << std::endl;
    }
    fout.close();

    fout.open(fnames[2]);
    for (unsigned i = 0; i < indices.size(); i++) {
        fout << indices.at(i) << " " << Times.at(i) << std::endl;
    }
    fout.close();

    return 0;
}
