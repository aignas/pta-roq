#include <iostream>
#include <vector>
#include <string>

#include "../libs/iocsv.hh"
#include "../libs/pulsar.hh"

int test_csvPulsar () {
    std::vector<Pulsar> pulsars, pulsars_read;
    std::string filename = "data_test/test_csv2pulsar.csv";

    const unsigned int pulsarNumber = 36;
    std::vector<double> range = {1e20, 1e24, 0, 3.1, -3.1, 3.1};

    // Randomize the pulsar structure and generate a schedule
    pulsarGrid::randomizeData (pulsars, pulsarNumber, range, 5e-9);

    pulsar2csv (filename, pulsars, ",");

    csv2pulsar (filename, pulsars_read, ",");

    int r = 0;

    // Check the pulsar size only at the moment
    if (pulsars_read.size() != 36) {
        // FIXME check whether the parameters of the two arrays match.
        r++;
    }

    return r;
}

int test_csvSources () {
    int r = 0;

    if (true) {
        r++;
    }

    return r;
}

int test_csvSchedule () {
    int r = 0;

    if (true) {
        r++;
    }

    return r;
}
