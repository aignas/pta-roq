
#include <vector>
#include <iostream>

#include "pulsar.hh"
#include "linalg.hh"
#include "random-helper.hh"
#include "signal-sampling.hh"

#include "reduced-basis.hh"

#include <omp.h>

#define CHUNKSIZE 20;

bool arrayContainsLong (long N, std::vector<long> & A);

void greedyReducedBasis (const unsigned long N,
                         void (*getData)(unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out),
                         std::vector<double> & A,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_param_out,
                         std::vector<std::vector<double> > & RB_out,
                         std::vector<double> & sigma_out) {
    // Initialise an empty Grammaian, the parameter space and calculate the size of the
    // parameter space (total number of points)

    std::vector<double> sigmaTrial (N, 0);
    std::vector<std::vector<double> > projectionCoeffs (N);
    
    // Seed choice (arbitrary). We just randomize the first choice. We also select the
    // first error arbitrary. This is just to have the same dimensions of two arrays
    std::vector<double> sigma_trial (N, 0);

    //std::vector<long> RB_N {random_uniform_int(0, N)};
    std::vector<long> RB_N {0};

    // Randomize the first basis
    sigma_out.resize(0);
    RB_out.resize(0);
    RB_param_out.resize(0);

    // Allocate some memory to params_trial and data arrays
    std::vector<double> params_trial, data;

    // Get data and parameters
    (*getData)(RB_N[0], params_trial, data);

    // Push back the values
    sigma_out.push_back(1000);
    RB_param_out.push_back(params_trial);
    RB_out.push_back(data);

    // Lets have a counter, which would mean that we search for basis even if epsilon
    // goes very small
    unsigned counter = 0;

    // The parameter space is large, so the computation will be expensive
    while (sigma_out.back() > epsilon and counter < 3) {
        // Construct the Gram matrix and its inverse
        // FIXME Use the previous Grammian to just extend it?
        std::vector<double> Grammian_inv;
        constructGrammian (Grammian_inv, RB_out, A);
        inverse (Grammian_inv);

        unsigned long chunk = CHUNKSIZE;

#pragma omp parallel for shared(sigma_trial, chunk) private(data, params_trial)
        // Stupidly traverse the entire parameter space
        // Can we edit the ranges where we are searching by discarding regions in
        // parameter space when we find the vectors? (Suggestion by Priscilla)
        for (unsigned long j = 0; j < N; j++) {
            sigma_trial.at(j) = 0;

            // Do not include the same point in the calculation
            if (arrayContainsLong(j, RB_N)) {
                continue;
            }

            (*getData)(j, params_trial, data);

            projectionResidual (data, RB_out, A, Grammian_inv, projectionCoeffs.at(j));

            sigma_trial.at(j) = innerProduct (data, A, data);
        }
        
        sigma_out.push_back(0);
        RB_N.push_back(0);

        findMax(sigma_trial, RB_N.back(), sigma_out.back());

        if (sigma_out.back() > epsilon) {
            // Add the lambda_i, which was found by maximizing the error
            (*getData)(RB_N.back(), params_trial, data);

            RB_param_out.push_back(params_trial);
            RB_out.push_back(data);
        } else {
            sigma_out.pop_back();
            RB_N.pop_back();
            counter++;
        }
    }

    // Free the memory occupied by the PRNG
    random_uniform_free();
}

void idToList (unsigned long idx, 
               std::vector<unsigned int> dim, 
               std::vector<unsigned int> & list_out) {
    // Clear the output vector
    list_out.resize(dim.size());

    // Construct the test vector of the parameters
    for (unsigned int i = 0; i < dim.size(); i++) {
        list_out[i] = idx % dim[i];
        idx = idx / dim[i];
    }
}

bool arrayContainsLong (long N, std::vector<long> & A) {
    bool r = false;
    for (unsigned i = 0; i < A.size(); i++) {
        if (A[i] == N) {
            r = true;
            break;
        }
    }
    
    return r;
}
