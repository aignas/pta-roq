#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

#include "pulsar.hh"

double random_simple (double LO, double HI) {
    return LO + (double)rand()/((double)RAND_MAX/(HI-LO));
}

void Pulsar::setRedNoise (double N, double f) {
    mRedNoise.resize(2);
    mRedNoise[0] = N;
    mRedNoise[1] = f;
}

void Pulsar::setPowerLawNoise (double N, double gamma) {
    mPowerLawNoise.resize(2);
    mPowerLawNoise[0] = N;
    mPowerLawNoise[1] = gamma;
}

namespace PulsarGrid {
    void randomizeData (std::vector<Pulsar> &Grid, unsigned int N, dvec range, double wnoise) {
        // Seed the random number generator
        srand (time(NULL));

        // Set the size of the containers
        Grid.resize(N);

        for (int i = 0; i < N; i++) {
            // Randomize angles in some range
            Grid[i].setAngles(random_simple(range[2], range[3]),  // random \theta
                              random_simple(range[4], range[5])); // random \phi

            // As well as the distance
            Grid[i].setDistance(random_simple(range[0], range[1]));

            // Randomize the noise values, such that all pulsars have different noises.
            Grid[i].setWhiteNoise(random_simple(0,wnoise));
        }
    }

    void generateSchedule (std::vector<Pulsar> &Grid, std::vector<double> initialTimes, double tFinal, double tMin, double tMax) {
        bool collectData = true;
        int N = initialTimes.size();
        std::vector<double> Times;

        // Copy the structure of the array for the time log
        std::vector<std::vector<double> > dates (N);

        // Iterate over each pulsar
        for (unsigned int i = 0; i < N; i++) {

            // Reset the time
            double t = initialTimes[i];

            // Generate a time series
            while (t < tFinal) {
                Times.push_back(t);

                // Advance in time
                t += random_simple(tMin, tMax);
            }

            Grid[i].setSchedule(Times);
        }
    }
}

