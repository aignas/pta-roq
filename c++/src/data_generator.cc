#include "linalg.hh"
#include "signal-model.hh"
#include "signal-sampling.hh"
#include "signal-vectors.hh"
#include "pulsar.hh"

#include <vector>
#include <iostream>
#include <fstream>

int main() {
    unsigned int pulsarNumber = 36;
    std::vector<double> range (6);
    range[0] = 50, range[0] = 100;
    range[1] = 0, range[0] = _M_PI;
    range[2] = -_M_PI, range[0] = _M_PI;

    // Define some variables for bigger timescales
    const double week = 7 * 3600 * 24,
                 year = 52 * week,
                 t_final = 15 * year,
                 dt_min = 2 * week,
                 dt_max = 2 * week;
    std::vector<double> t_init (pulsarNumber, 0),
                        Times;
    std::vector<unsigned short> indices;

    std::vector<Pulsar> pulsars;

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Randomizing the initial conditions" << std::endl;
    pulsarGrid::randomizeData (pulsars, pulsarNumber, range, 0.00005);
    std::cout << "Randomizing the initial conditions" << std::endl;
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, indices, Times);

    // Add a single source
    std::vector<std::vector<double> > sources (1);
    sources[0].resize(8);
    sources[0][0] = 1;
    sources[0][1] = 1e31;
    sources[0][2] = 0;
    sources[0][3] = 0;
    sources[0][4] = 0;
    sources[0][5] = 0.3;
    sources[0][6] = 0.5;
    sources[0][7] = 1e-8;

    // Generate a residual vector
    std::cout << "Generating a residual vector" << std::endl;
    std::vector<double> r;
    generateSample (r, pulsars, indices, Times, sources);

    // Output the stuff
    // FIXME

    return 0;
}
