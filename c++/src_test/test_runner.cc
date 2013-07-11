#include <iostream>
#include <string>

#include "linalg.hh"
#include "reduced-basis.hh"

void testing_wrap( int (*f)(void), std::string name, unsigned int & r, unsigned int & t);
void title_wrap (std::string name);

int main () {
    std::cout << "Launching the testing suite" << std::endl;

    unsigned int r = 0, t = 0;

    title_wrap("linalg");
    testing_wrap(test_arrayEqual, "arrayEqual", r, t);
    testing_wrap(test_dotProduct, "dotProduct", r, t);
    testing_wrap(test_innerProduct, "innerProduct", r, t);
    testing_wrap(test_inverse, "inverse", r, t);
    testing_wrap(test_constructGrammian, "constructGrammian", r, t);
    testing_wrap(test_projectionResidual, "projectionResidual", r, t);
    testing_wrap(test_linspace, "linspace", r, t);
    testing_wrap(test_findMax, "findMax", r, t);

    title_wrap("reduced-basis");
    testing_wrap(test_idToList, "idToList", r, t);
    testing_wrap(test_greedyReducedBasis, "greedyReducedBasis", r, t);

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
