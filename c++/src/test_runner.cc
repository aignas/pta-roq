#include <iostream>
#include <string>

#include "test/signal/model.hh"
#include "test/linalg.hh"
#include "test/roq.hh"
#include "test/iocsv.hh"

void testing_wrap( int (*f)(void), std::string name, unsigned int & r, unsigned int & t);
void title_wrap (std::string name);

int main () {
    std::cout << "\nLaunching the testing suite" << std::endl;

    unsigned int r = 0, t = 0;

    title_wrap("signal/model");
    // This test is not right, because of floating point number addition in the formula
    // pulsarTerm. That formula shouldn't be used;
    //testing_wrap(test_pulsarTerm, "pulsarTerm", r, t);

    title_wrap("linalg");
    testing_wrap(test_arrayEqual, "arrayEqual", r, t);
    testing_wrap(test_dotProduct, "dotProduct", r, t);
    testing_wrap(test_innerProduct, "innerProduct", r, t);
    testing_wrap(test_inverse, "inverse", r, t);
    testing_wrap(test_constructGrammian, "constructGrammian", r, t);
    testing_wrap(test_projectionResidual, "projectionResidual", r, t);
    testing_wrap(test_linspace, "linspace", r, t);
    testing_wrap(test_findMax, "findMax", r, t);
    testing_wrap(test_findMin, "findMin", r, t);

    title_wrap("reduced-basis");
    testing_wrap(test_idToList, "idToList", r, t);
    testing_wrap(test_greedyReducedBasis, "greedyReducedBasis", r, t);
    std::cout << "The following tests are not implemented yet!" << std::endl;
    testing_wrap(test_greedyEIMpoints, "greedyEIMpoints", r, t);
    testing_wrap(test_constructROQ, "constructROQ", r, t);

    title_wrap("iocsv");
    testing_wrap(test_csvPulsar,   "csv2pulsar and pulsar2csv",     r, t);
    testing_wrap(test_csvSources,  "csv2sources and sources2csv",   r, t);
    testing_wrap(test_csvSchedule, "csv2schedule and schedule2csv", r, t);

    std::cout << "\nTest results:\n\t" << r << " out of " << t << " tests were failed" << std::endl;

    return 0;
}

void testing_wrap( int (*f)(void), std::string name, unsigned int & r, unsigned int & t) {
    // Increment the number of tests run;
    t++;

    // Increment the number of tests failed
    int k = (*f)();
    r += k;

    std::string markup, result;


    if (k==1) {
        result = "failed";
        markup = "\x1b[31;1m";
    } else {
        result = "passed";
        markup = "\x1b[32;0m";
    }

    std::cout << markup << "Testing " << name << ": " << result << "\x1b[0m" << std::endl;
}

void title_wrap(std::string name) {
    std::cout << "\x1b[1m" << "\nTesting the " << name << ".hh functions" << "\x1b[0m" << std::endl;
}
