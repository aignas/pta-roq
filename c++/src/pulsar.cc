#include <vector>
#include <pulsar.hh>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

random (double LO, double HI) {
    return LO + (double)rand()/((double)RAND_MAX/(HI-LO));
}

PulsarGrid::PulsarGrid (unsigned int N, dvec range, double wnoise) {
    PulsarGrid::randomizeData(N, range, noise);
}

PulsarGrid::randomizeData (unsigned int N, dvec range, double wnoise) {
    // Seed the random number generator
    srand (time(NULL));

    // Set the size of the containers
    PulsarGrid::setSize(N);

    for (int i = 0; i < N; i++) {
        // Randomize angles in some range
        mAngle[i][0] = random(range[2], range[3]);
        mAngle[i][1] = random(range[4], range[5]);

        // As well as the distance
        mDistance[i] = random(range[0], range[1]);

        // Randomize the noise values, such that all pulsars have different noises.
        mWhiteNoise[i] = random(0,wnoise);
    }
}

PulsarGrid::setAngles (std::vector<dvec> angles) {
    mAngle = angles;
}

PulsarGrid::setSize (unsigned int N) {
    mAngle.resize(N);
    mDistance.resize(N);
    mWhiteNoise.resize(N);

    for (int i = 0; i < N; i++) {
        // Randomize angles in some range
        mAngle[i].resize(2);
    }
}

PulsarGrid::setWhiteNoise (dvec wnoise) {
    mWhiteNoise = wnoise;
}
