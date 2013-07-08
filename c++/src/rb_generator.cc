
#include <vector>
#include <valarray>
#include <iostream>
#include <fstream>

#include "pulsar.hh"
#include "linalg.hh"
#include "signal-sampling.hh"

#include "reduced-basis.hh"

int main() {
    unsigned int pulsarNumber = 6;
    std::vector<double> range (6);
    range[0] = 50, range[0] = 100;
    range[1] = 0, range[0] = _M_PI;
    range[2] = -_M_PI, range[0] = _M_PI;

    // Define some variables for bigger timescales
    const double week = 7 * 3600 * 24,
                 year = 52 * week,
                 t_final = 0.5 * year,
                 dt_min = 2 * week,
                 dt_max = 2 * week;

    std::vector<double> t_init (pulsarNumber, 0),
                        Times;
    std::vector<unsigned short> indices;

    std::vector<Pulsar> pulsars;

    // Randomize the pulsar structure and generate a schedule
    std::cout << "Randomizing the initial conditions" << std::endl;
    pulsarGrid::randomizeData (pulsars, pulsarNumber, range, 0.00005);
    std::cout << "Constructing the schedule of measuremnts" << std::endl;
    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, indices, Times);

    // Construct the parameter space
    std::vector<std::vector<double> > params (8);
    linspace(params[0], 1,1,1);
    linspace(params[1], 1e9, 1e9, 1);
    linspace(params[2], 0.1, 0.1, 1);
    linspace(params[3], 0.1, 0.1, 1);
    linspace(params[4], 0.1, 0.1, 1);
    linspace(params[5], 0.5, 1.0, 10);
    linspace(params[6], 0.0, 1.0, 10);
    linspace(params[7], 9.9e-9, 1.1e-8, 10);

    // Set the error
    const double epsilon = 1e-4;

    // Constructing a covariance matrix
    std::cout << "Constructing the covariance matrix" << std::endl;
    std::vector<double> C;
    genCovarianceMatrix (C, pulsars, indices, Times, true, false, false, false);

    // Generate a residual vector
    std::cout << "Generating the reduced basis set" << std::endl;
    std::vector<std::vector<double> > RB;
    greedyReducedBasis (indices, Times, pulsars,
            params, C, epsilon, RB);

    return 0;
}
