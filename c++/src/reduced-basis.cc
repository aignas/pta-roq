
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
                         std::vector<double> & sigma_out,
                         bool verbose) {
    // Initialise an empty Grammaian, the parameter space and calculate the size of the
    // parameter space (total number of points)

    std::vector<double> sigmaTrial (N, 0),
                        templateNorm (N, 0);
    std::vector<std::vector<double> > RB_out_hat;
    
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
    RB_out_hat.push_back(data);
    matrixVectorProduct (A, RB_out.back(), RB_out_hat.back());

    // Output if demanded
    if (verbose) {
        std::cout << "N, id, error, p1, p2, p3, p4, p5, p6, p7, p8" << std::endl;
        std::cout << sigma_out.size() << ", " << RB_N.back() << ", " << sigma_out.back();
        for (unsigned j = 0; j < RB_param_out.back().size(); j++) {
            std::cout << ", " << RB_param_out.back().at(j);
        }
        std::cout << std::endl;
    }
    // Lets have a counter, which would mean that we search for basis even if epsilon
    // goes very small
    unsigned counter = 0;

    std::vector<double> Grammian, Grammian_inv;

    // Calculate the template norms.
    // FIXME This will need to be stored in a file for likelyhood evaluations,
    // because in this way I can ommit expensive inner products there and here
    // I can't get away without calculating them, so it's a win/win situation.
#pragma omp parallel for shared(sigma_trial, chunk) private(data, params_trial)
    for (unsigned long i = 0; i < N; i++) {
        (*getData)(i, params_trial, data);
        templateNorm.at(i) = innerProduct (data, A, data);
    }

    // The parameter space is large, so the computation will be expensive
    while (sigma_out.back() > epsilon) {
        // Construct the Gram matrix and its inverse
        // FIXME Use the previous Grammian to just extend it?
        extendGrammianOptimized (Grammian, Grammian_inv, RB_out, RB_out_hat);
        //constructGrammianOptimized (Grammian, RB_out, RB_out_hat);

        unsigned long chunk = CHUNKSIZE;

        // FIXME Stupidly traverse the entire parameter space. Possible
        // improvements:
        //  * Do an MCMC search
        //
        // parallelize it with OpenMP
#pragma omp parallel for shared(sigma_trial, templateNorm, chunk) private(data, params_trial)
        for (unsigned long j = 0; j < N; j++) {
            // Reset the error to zero
            sigma_trial.at(j) = 0;

            // Do not include the same point in the calculation
            if (arrayContainsLong(j, RB_N)) {
                continue;
            }

            (*getData)(j, params_trial, data);

            
            sigma_trial.at(j) = projectionErrorStable (templateNorm[j], data, 
                    RB_out_hat, Grammian, Grammian_inv);
        }
        
        sigma_out.push_back(0);
        RB_N.push_back(0);

        findMax(sigma_trial, RB_N.back(), sigma_out.back());

        // Add the lambda_i, which was found by maximizing the error
        (*getData)(RB_N.back(), params_trial, data);

        RB_param_out.push_back(params_trial);
        RB_out.push_back(data);
        RB_out_hat.push_back(data);
        // Multiply one of the basis vectors with the covariance matrix
        matrixVectorProduct (A, RB_out.back(), RB_out_hat.back());

        // Output if demanded
        if (verbose) {
            std::cout << sigma_out.size() << ", " << RB_N.back() << ", " << sigma_out.back();
            for (unsigned j = 0; j < RB_param_out.back().size(); j++) {
                std::cout << ", " << RB_param_out.back().at(j);
            }
            std::cout << std::endl;
        }
    }

    // Free the memory occupied by the PRNG
    random_uniform_free();
}

void idToList (unsigned long idx, 
               std::vector<unsigned int> & dim, 
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
