
#include <vector>
#include <valarray>
#include <iostream>

#include "pulsar.hh"
#include "linalg.hh"
#include "random-helper.hh"
#include "signal-sampling.hh"

#include "reduced-basis.hh"

void greedyReducedBasis (std::vector<unsigned short> & indices,
                         std::vector<double> & schedule,
                         std::vector<Pulsar> & pulsars,
                         std::vector< std::vector<double> > & params,
                         std::vector<double> & covarianceMatrix,
                         const double & epsilon,
                         std::vector<std::vector<double> > & RB_out) {
    // Initialise an empty Grammaian, the parameter space and calculate the size of the
    // parameter space (total number of points)

    unsigned long totalNumber = 1;
    for (unsigned int i = 0; i < params.size(); i++) {
        totalNumber *= params.at(i).size();
    }

    std::vector<double> sigmaTrial (totalNumber, 0);
    std::vector<std::vector<double> > projectionCoeffs (totalNumber);
    
    // Seed choice (arbitrary). We just randomize the first choice. We also select the
    // first error arbitrary. This is just to have the same dimensions of two arrays
    std::vector<double> sigma (1, 1);
    std::vector<double> sigma_trial (totalNumber, 0);

    // Randomize the first basis
    RB_out.resize(0);
    // Allocate some memory to params_trial
    std::vector<std::vector<double> > params_trial;
    std::cout << random_uniform_int(0, totalNumber, 1) << std::endl;

    idToParam(random_uniform_int(0, totalNumber), params, params_trial);

    RB_out.push_back(params_trial[0]);

    std::vector<double> data (indices.size(), 0);

    std::vector<std::vector<double> > RB_Training;
    generateSample (data, pulsars, indices, schedule, RB_out);
    RB_Training.push_back(data);

    // The parameter space is large, so the computation will be expensive
    while (sigma.back() > epsilon) {
        // Construct the Gram matrix and its inverse
        // FIXME 
        std::vector<double> Grammian;
        constructGrammian (Grammian, RB_Training, covarianceMatrix);
        std::vector<double> Grammian_inv = Grammian;
        inverse (Grammian_inv);

        // FIXME Make this parallel with the pragmas
        // Stupidly traverse the entire parameter space
        // NOTE: We could have MCMC or a simple MC method as well
        // Can we edit the ranges where we are searching by discarding regions in
        // parameter space when we find the vectors? (Suggestion by Priscilla)
        for (unsigned long i = 0; i < totalNumber; i++) {

            // the params_trial does not have to be zeroed away
            // This is a temporary variable
            std::vector<std::vector<double> > params_tmp;
            idToParam(i, params, params_tmp);

            // Calculate the data
            generateSample (data, pulsars, indices, schedule, params_tmp);

            projectionResidual (data, RB_Training, 
                                covarianceMatrix, Grammian_inv, 
                                projectionCoeffs.at(i));

            sigma_trial.at(i) = norm (data, covarianceMatrix);

            /*
            for (unsigned j = 0; j < params_tmp[0].size(); j++) {
                std::cout << params_tmp[0][j] << "\t";
            }
            std::cout << sigma_trial.at(i) << std::endl;
            */
        }

        // Find a maximum in the calculated errors
        double sigma_max = sigma_trial.at(0),
               sigma_max_tmp;
        unsigned long sigma_arg_max = 0;
        for (unsigned long i = 1; i < totalNumber; i++) {
            sigma_max_tmp = sigma_trial.at(i);
            if (sigma_max < sigma_max_tmp) {
                sigma_max = sigma_max_tmp;
                sigma_arg_max = i;
            }
        }

        sigma.push_back(sigma_max);

        // Add the lambda_i, which was found by maximizing the error
        idToParam(sigma_arg_max, params, params_trial);
        RB_out.push_back(params_trial[0]);

        generateSample (data, pulsars, indices, schedule, params_trial);
        RB_Training.push_back(data);

        std::cout << "Number of RB and the error: "
                  << sigma.size() << " " << sigma.back() << std::endl;
    }

    // Free the memory
    random_uniform_free();
}

void idToParam (unsigned long idx, 
                std::vector<std::vector<double> > & space, 
                std::vector<std::vector<double> > & param_out) {

    std::vector<double> tmp (space.size(), 0);
    param_out.clear();

    // Construct the test vector of the parameters
    for (unsigned int i = 0; i < tmp.size(); i++) {
        unsigned int n = idx % space.at(i).size();
        tmp.at(i) = space.at(i).at(n);
        idx = idx / space[i].size();
    }

    // This could be generalized to search of more than one source
    param_out.push_back(tmp);
}

