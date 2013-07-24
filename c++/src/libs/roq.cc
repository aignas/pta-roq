#include <vector>
#include <iostream>
#include <omp.h>

#include "linalg.hh"
#include "random-helper.hh"

#include "roq.hh"

#define CHUNKSIZE 10;

bool arrayContainsLong (long N, std::vector<long> & A);

void greedyReducedBasis (const unsigned long N,
                         void (*getData)(unsigned long idx, std::vector<double> & params_out, std::vector<double> & data_out),
                         std::vector<double> & A,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_param_out,
                         std::vector<std::vector<double> > & RB_out,
                         std::vector<double> & G_out,
                         std::vector<double> & templateNorm_out,
                         std::vector<double> & sigma_out,
                         bool verbose) {
    // Initialise an empty Grammaian, the parameter space and calculate the size of the
    // parameter space (total number of points)

    std::vector<double> sigmaTrial (N, 0);
    std::vector<std::vector<double> > RB_out_hat,
                                      coefficientsTmp (N);
    
    // Seed choice (arbitrary). We just randomize the first choice. We also select the
    // first error arbitrary. This is just to have the same dimensions of two arrays
    std::vector<long> RB_N {random_uniform_int(0, N)};
    //std::vector<long> RB_N {0};

    // Resize the output containers
    sigma_out.clear();
    RB_out.clear();
    RB_param_out.clear();
    G_out.clear();
    templateNorm_out.resize(N);

    if (verbose) {
        std::cout << "Generating the reduced basis: " << std::endl;
    }
    // Allocate some memory to paramsTrial and data arrays
    std::vector<double> paramsTrial, data;

    // Get data and parameters
    (*getData)(RB_N[0], paramsTrial, data);

    // Push back the values
    sigma_out.push_back(100);
    RB_param_out.push_back(paramsTrial);
    RB_out.push_back(data);
    RB_out_hat.push_back(data);
    matrixVectorProduct (A, RB_out.back(), RB_out_hat.back());

    // Calculate the template norms.
    // FIXME This will need to be stored in a file for likelyhood evaluations,
    // because in this way I can ommit expensive inner products there and here
    // I can't get away without calculating them, so it's a win/win situation.
    if (verbose) {
        std::cout << "Calculating the norms of the signal templates: "; std::cout.flush();
    }

    unsigned int chunk = CHUNKSIZE;

#pragma omp parallel shared(templateNorm_out, chunk) private(data, paramsTrial)
    {
#pragma omp for
    for (unsigned long i = 0; i < N; i++) {
        (*getData)(i, paramsTrial, data);
        templateNorm_out.at(i) = innerProduct (data, A, data);
    }
    } 

    std::vector<double> Grammian_inv;
    if (verbose) {
        std::cout << "DONE" << std::endl;
        std::cout << "The relative error at each iteration:" << std::endl;
    }

    // The parameter space is large, so the computation will be expensive
    while (sigma_out.back() > epsilon) {
        // // Normalize the basis vectors:
        // double alpha = sqrt(templateNorm.at(RB_N.back()));
        // for (unsigned int j = 0; j < RB_out.back().size(); j++) {
        //     RB_out.back()[j] /= alpha;
        //     RB_out_hat.back()[j] /= alpha;
        // }

        // Construct the Gram matrix and its inverse
        // FIXME Use the previous Grammian to just extend it?
        try {
            extendGrammianOptimized (G_out, RB_out, RB_out_hat);
        } catch (const char *msg) {
            std::cerr << msg << std::endl;
        }

        // Check if the grammian is invertible
        int i = inverseATLAS(G_out, Grammian_inv);
        if (i != 0) {
            // This part is questionable, but let's see, how it goes.
            long argdel;
            double sigma_min = 0;
            sigma_out.erase(sigma_out.begin());
            findMin(sigma_out, argdel, sigma_min);

            argdel++;

            sigma_out.erase(sigma_out.begin() + argdel, sigma_out.end());
            RB_N.erase(RB_N.begin() + argdel, RB_N.end());
            RB_out.erase(RB_out.begin() + argdel, RB_out.end());
            RB_param_out.erase(RB_param_out.begin() + argdel, RB_param_out.end());

            break;
        }

        // Precalculate the projections of the signal templates onto the basis vectors
        // (do not worry about the orthogonality), since these should stay constant, we
        // can reuse them
#pragma omp parallel shared(RB_out_hat, coefficientsTmp, chunk) private(data, paramsTrial)
        {
#pragma omp for
            for (unsigned long j = 0; j < N; j++) {
                (*getData)(j, paramsTrial, data);

                coefficientsTmp.at(j).push_back(dotProduct(data, RB_out_hat.back()));
            }
        }

        // FIXME Stupidly traverse the entire parameter space. Possible
        // improvements:
        //  * Do an MCMC search
        //
        // parallelize it with OpenMP
#pragma omp parallel shared(sigmaTrial, templateNorm_out, chunk, G_out, Grammian_inv, RB_out_hat, RB_N) private(data, paramsTrial)
        {
#pragma omp for
            for (unsigned long j = 0; j < N; j++) {
                // Reset the error to zero
                sigmaTrial.at(j) = 0;

                (*getData)(j, paramsTrial, data);


                sigmaTrial.at(j) = projectionErrorStableNew (templateNorm_out.at(j), coefficientsTmp.at(j), Grammian_inv);
            }
        }

        // Add the max and argmax of sigmaTrial
        sigma_out.push_back(0);
        RB_N.push_back(0);
        findMax(sigmaTrial, RB_N.back(), sigma_out.back());

        if (verbose) {
            std::cout << "\t" << sigma_out.size() << "\t" << sigma_out.back() << std::endl;
        }

        // Add the lambda_i, which was found by maximizing the error
        (*getData)(RB_N.back(), paramsTrial, data);

        RB_param_out.push_back(paramsTrial); RB_out.push_back(data); 
        // Multiply one of the basis vectors with the covariance matrix
        RB_out_hat.push_back(data);
        matrixVectorProduct (A, data, RB_out_hat.back());
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

void reduceOrder (std::vector<long> & idx,
                  std::vector<double> & h,
                  std::vector<double> & h_out) {
    unsigned int N = idx.size();
    h_out.clear();

    // Construct \vec{h}
    for (unsigned i = 0; i < N; i++) {
        h_out.push_back(h.at(idx.at(i)));
    }
}

void constructInterpolant (std::vector<std::vector<double> > & RB,
                           std::vector<double> & h,
                           std::vector<long> & idx,
                           std::vector<double> & A_out,
                           std::vector<double> & data_out) {
    // One of the matrices used here:
    unsigned int N = idx.size();
    std::vector<double> h_tmp(N), c_tmp(N);

    A_out.clear();
    A_out.resize(N*N);

    // Construct the matrix
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            A_out.at(i+j*N) = RB.at(i).at(idx.at(j));
        }
    }

    // Invert A
    inverseATLASOverwrite(A_out);

    reduceOrder (idx, h, h_tmp);

    // Calculate \vec{c}
    matrixVectorProduct(A_out,h_tmp,c_tmp);

    // Calculate the last product:
    data_out.clear();
    data_out.resize(h.size());

    for (unsigned i = 0; i < N; i++) {
        axpyProduct(c_tmp.at(i), RB.at(i), data_out);
    }
}

void greedyEIMpoints (std::vector<std::vector<double> > & RB_param,
                      std::vector<std::vector<double> > & RB,
                      std::vector<long> & idx_out,
                      std::vector<double> & points_out,
                      std::vector<double> & A_out) {
    unsigned long N = RB.size(), Ndat = RB.at(0).size();
    std::vector<double> error (Ndat);

    // Clear the output arrays:
    idx_out.clear();
    points_out.clear();

    // Initiate the seed value
    idx_out.push_back(0);
    points_out.push_back(0);
    error = RB.at(0);

    // Find a modulus
    for (unsigned i=0; i < error.size(); i++) {
        error[i] = fabs(error[i]);
    }

    findMax(error, idx_out.back(), points_out.back());

    for (unsigned int i = 1; i < N; i++) {
        // Construct the interpolant
        constructInterpolant (RB, RB.at(i), idx_out, A_out, error);
        // Find the error
        axpyProduct(-1, RB.at(i), error);

        idx_out.push_back(0);
        points_out.push_back(0);

        // Find a modulus
        for (unsigned j=0; j < error.size(); j++) {
            error[j] = fabs(error[j]);
        }

        findMax(error, idx_out.back(), points_out.back());
    }
}

void constructROQ (std::vector<double> & r,
                   std::vector<double> & r_out,
                   std::vector<long> & idx,
                   std::vector<double> & G,
                   std::vector<double> & A) {
    unsigned int N = idx.size();

    // One of the matrices used here:
    std::vector<double> r_tmp (N, 0),
                        r_out_tmp (N,0);
    r_out.resize(N);

    // do the products
    reduceOrder (idx, r, r_tmp);
    matrixVectorProduct (G, r_tmp, r_out_tmp);

    matrixTranspVectorProduct (A, r_out_tmp, r_out);
    std::cout << "Debug" << std::endl;
}

