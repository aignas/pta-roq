#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "output-helper.hh"

void outputVector (std::vector<double> & A, std::string separator) {
    const unsigned int N = A.size();

    std::cout << A.at(0);
    for (unsigned int i = 1; i < N; i++) {
        std::cout << separator << A.at(i);
    }
}

void outputVector (std::vector<int> & A, std::string separator) {
    const unsigned int N = A.size();

    std::cout << A.at(0);
    for (unsigned int i = 1; i < N; i++) {
        std::cout << separator << A.at(i);
    }
}

void outputMatrix (std::vector<double> & A, std::string sep1, std::string sep2) {
    const unsigned int N = sqrt(A.size());

    for (unsigned int i = 0; i < N; i++) {
        std::cout << A.at(i*N);
        for (unsigned int j = 1; j < N; j++) {
            std::cout << sep1 << A.at(j + N*i);
        }
        if (i != N-1) {
            std::cout << sep2;
        }
    }
}

void outputMatrix (std::vector<int> & A, std::string sep1, std::string sep2) {
    const unsigned int N = sqrt(A.size());

    for (unsigned int i = 0; i < N; i++) {
        std::cout << A.at(i*N);
        for (unsigned int j = 1; j < N; j++) {
            std::cout << sep1 << A.at(j + N*i);
        }
        if (i != N-1) {
            std::cout << sep2;
        }
    }
}
