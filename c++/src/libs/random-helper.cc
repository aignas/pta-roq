
#include "gsl/gsl_rng.h"
#include "random-helper.hh"

// This is a simple PRNG with a range of [0;1] using the gsl library
double random_uniform (const double LO, 
                       const double HI, 
                       const unsigned long int seed) {
    // initialize a pointer rng
    static gsl_rng* rng = 0;

    if ( rng == 0 ) {
        // haven't used the function yet
        rng = gsl_rng_alloc( gsl_rng_default );
    }

    if ( seed != 0 ) {
        // need to set the seed
        gsl_rng_set( rng, seed );
    }

    // return the prn
    return LO + (double) gsl_rng_uniform(rng)*(HI-LO);
}

// This is a simple PRNG with a range of [0;1] using the gsl library
long random_uniform_int (const long LO, 
                         const long HI, 
                         const unsigned long int seed) {
    // initialize a pointer rng
    static gsl_rng* rng = 0;

    if ( rng == 0 ) {
        // haven't used the function yet
        rng = gsl_rng_alloc( gsl_rng_default );
    }

    if ( seed != 0 ) {
        // need to set the seed
        gsl_rng_set( rng, seed );
    }

    // return the prn
    return LO + (long) gsl_rng_uniform_int(rng, HI-LO);
}

void random_uniform_free () {
    // initialize a pointer rng
    static gsl_rng* rng = 0;

    if ( rng != 0 ) {
        gsl_rng_free(rng);
    }
}

