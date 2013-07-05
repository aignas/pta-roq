/**
 * The random number generator helper library
 */

/**
 * Generate a random number between [0,1] if the limits are unspecified
 *
 * @param LO The low value of the range (default 0)
 * @param HI The high value of the range (default 1)
 * @param seed The seed value (default 0)
 */
double random_uniform (const double LO = 0, 
                       const double HI = 1, 
                       const unsigned long int seed = 0);

/**
 * Free the random number generator
 */
void random_uniform_free (void);
