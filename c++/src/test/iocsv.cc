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
    if (pulsars_read.size() != pulsarNumber) {
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
    // We need some pulsar data
    test_csvPulsar ();

    std::vector<Pulsar> pulsars;
    std::string filename_pulsar = "data_test/test_csv2pulsar.csv",
                filename_sched = "data_test/test_csv2schedule.csv";

    // Read some pulsar data
    csv2pulsar (filename_pulsar, pulsars, ",");

    double t_final = 5*52, dt_min = 2, dt_max = 2;
    std::vector<unsigned short> indices, indices_read;

    const unsigned int pulsarNumber = 36;
    std::vector<double> Times, Times_read, t_init (pulsarNumber, 0);

    pulsarGrid::generateSchedule (t_init, t_final, dt_min, dt_max, indices, Times);

    arraysShortDouble2csv (filename_sched, indices, Times, ",");
    csv2arraysShortDouble (filename_sched, indices_read, Times_read, ",");

    int r = 0;

    if (indices.size() != indices.size() or Times.size() != Times_read.size()) {
        std::cerr << indices.size() << " != " << indices.size() 
                  << " or " << Times.size() << " != " << Times_read.size()
                  << std::endl;
        r++;
    }

    return r;
}
