#include <vector>
#include <iostream>
#include <cmath>
#include <ctime>

#include "random-helper.hh"
#include "linalg.hh"
#include "mcmc.hh"

/*=========== MAIN ===========*/
    // THIS MIGHT BE USEFUL (this was used for creating a directory in Dr. Gair's code)
    //system(shell_command);

    // initialization of run Parameters

    // Calc C^{-1}

    // Remmember!
    // search_start=clock();
    // Do MCMC
    // end=clock();

    // Data save and output statistics

double ComputeLogLikelihood (std::vector<double> & data,
                             std::vector<double> & signalTmp,
                             std::vector<double> & CovM) {
    return dotProduct(data,signalTmp) - norm(signalTmp, CovM)/2;
}

void drawFromPrior(std::vector<double> & parmin, 
                   std::vector<double> & parmax,
                   std::vector<double> & thprop) {
    // There needs to be a way to introduce different priors here.
    thprop.resize(parmin.size());

    // For each parameter space dimension generate a random number, i.e. select a point
    // randomly in the multi-D space
    for (unsigned int i = 0; i < parmin.size(); i++) {
        thprop[i] = random_uniform(parmin.at(i), parmax.at(i));
    }
}

void computeFisherMatrix(std::vector<double> & thnow,
                         std::vector<double> & fMatrix) {
    /*
     * Use the derivatives to calculate the fisher matrix. This might be however very
     * expensive as I need to have m^2 inner products, where m = number of parameters in
     * search and all of them scale as n^3, where n is the number of data points in the
     * vectors, which might be very large. This is not a good thing.
     */
}

void GaussianProposal (std::vector<double> & parmin, 
                       std::vector<double> & parmax,
                       std::vector<double> & thnow,
                       std::vector<double> & thprop) {
    drawFromPrior(parmin, parmax, thprop);
}

void mcmcProposal(std::vector<double> & parmin,
                  std::vector<double> & parmax, 
                  std::vector<double> & thnow, 
                  std::vector<double> & thprop,
                  bool & useGaussian) {
    if (useGaussian) {
        GaussianProposal(thnow, thprop, parmin, parmax);
    } else {
        drawFromPrior(parmin, parmax, thprop);
    }
}

void mcmcSearch(std::vector<double> & data,
                std::vector<double> & CovM,
                void (*getTmp) (std::vector<double> & thnow,
                                std::vector<double> & signalTmp),
                std::vector<std::vector<double> > & paramSpace, 
                const unsigned int NMCPoints,
                const bool dataSave,
                bool & ROQ) {
    double ratio, likenew, likeold;
    bool accept;
    // Various containers
    std::vector<double> parmin,
                        parmax,
                        thnow,
                        thprop,
                        signalTmp;

    std::clock_t start;
    start = std::clock();
    computeLogLikelihood(data, signalTmp, CovM);
    std::cout << "Likelihood is evaluated in " 
              << ( std::clock() - start ) / (double) CLOCKS_PER_SEC
              << " seconds." << std::endl;

    // Data saving

    // FIXME We can parallelize this in order to use more cores. One way would be to
    // parallelize the algo, another, just search the parameter space with a few chains
    // at the same time. They chains could explore the same or different parameter
    // spaces.

    // Metropolis-Hastings algo
    for (unsigned int i = 0; i < NMCPoints; i++) {
        mcmcProposal(parmin, parmax, thnow, thprop);

        getTmp(thprop, signalTmp);

        likenew = computeLogLikelihood(data, signalTmp, CovM);
        ratio = exp(likenew-likeold);


        if ((ratio > 1.) or (random_uniform(0,1) < ratio)) {
            likeold=likenew;
            axpyProduct(1, thprop, 0, thnow);
        }

        if (dataSave) {
            // Output to files. I should pass the file names as well.
            // Stuff to output is:
            //      - thnow,
            //      - likeold
        }
    }

    random_uniform_free();
}
