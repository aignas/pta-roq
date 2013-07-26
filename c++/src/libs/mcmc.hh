#include <vector>

double computeLogLikelihood (std::vector<double> & data,
                             std::vector<double> & signalTmp,
                             std::vector<double> & CovM);

/**
 * Draw from prior
 */
void drawFromPrior(std::vector<double> parmin, 
                   std::vector<double> parmax,
                   std::vector<double> & prop);

/**
 * Proposal wrapper function
 */
void mcmcProposal(std::vector<double> & thnow, 
                  std::vector<double> & thprop,
                  std::vector<double> & parmin,
                  std::vector<double> & parmax, 
                  bool & useGaussian = false);

/**
 * Search the parameter space with MCMC
 *
 * @param data The residual Can be either roq or non ROQ, the main structure of the problem
 *        remains the same. The residual must be premultiplied by an appropriate
 *        covariance matrix, which in ROQ case is not the same as CovM, but in the
 *        simple case it is the same matrix as CovM.
 * @param CovM The covariance matrix used to calculate the norm of the signal template.
 * @param getTmp This function is used to get transform the data of thnow and return the
 *        signal template. Since it is defined externally, we can make the mcmc
 *        independent of actual model used.
 * @param paramSpace The vector of vectors having all the important parameter line
 *        spaces in different rows
 * @param NMCPoints this is the total number of Monte Carlo points to be used in this
 *        simulation. This means, that after some number of trials, the Markov chain
 *        will be terminated.
 * @param Whether to save data into a file
 * @param ROQ True or false depending whether we are using the ROQ method or not.
 */
void mcmcSearch(std::vector<double> & data,
                std::vector<double> & CovM,
                void (*getTmp) (std::vector<double> & thnow,
                                std::vector<double> & signalTmp),
                std::vector<std::vector<double> > & paramSpace, 
                const unsigned int NMCPoints,
                const bool dataSave,
                const bool ROQ);
