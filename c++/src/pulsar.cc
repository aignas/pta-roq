#include <vector>
#include <cmath>

#include "pulsar.hh"
#include "random-helper.hh"

Pulsar::Pulsar () {
        mDistance = 0;
        mTheta = 0;
        mPhi = 0;
        mWhiteNoise = 0;
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

std::vector<double> Pulsar::getUnitVector () {
    std::vector<double> v (3, 0);

    v[0] = sin(mTheta);
    v[1] = v[0];
    v[2] = cos(mTheta);
    v[0] *= cos(mPhi);
    v[1] *= sin(mPhi);

    return v;
}

namespace pulsarGrid {
    void randomizeData (std::vector<Pulsar>& Grid, unsigned int N, dvec range, double wnoise) {
        for (unsigned int i = 0; i < N; i++) {
            Pulsar tmp;

            // Randomize angles in some range
            tmp.setAngles(random_uniform(range[2], range[3]),  // random \theta
                          random_uniform(range[4], range[5])); // random \phi

            // As well as the distance
            tmp.setDistance(random_uniform(range[0], range[1]));

            // Randomize the noise values, such that all pulsars have different noises.
            tmp.setWhiteNoise(random_uniform(0,wnoise));

            Grid.push_back(tmp);
        }

        random_uniform_free();
    }

    void generateSchedule ( std::vector<double>& initialTimes, 
                           double tFinal, double tMin, double tMax,
                           std::vector<unsigned short>& indices,
                           std::vector<double>& Times ) {
        unsigned int N = initialTimes.size();

        // Copy the structure of the array for the time log
        std::vector<std::vector<double> > dates (N);

        // Iterate over each pulsar
        for (unsigned int i = 0; i < N; i++) {

            // Reset the time
            double t = initialTimes.at(i);

            // Generate a time series
            while (t < tFinal) {
                Times.push_back(t);
                indices.push_back(i);

                // Advance in time
                t += random_uniform(tMin, tMax);
            }
        }

        random_uniform_free();
    }
}

