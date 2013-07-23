#include <vector>

// FIXME: Do I need this? because at the moment the pulsar grid is not used in this file
// and might never be used
#include "../pulsar.hh"

#ifndef _SIGNAL_MODEL_HXX_
#define _SIGNAL_MODEL_HXX_
// Define PI to a high precision
const double _M_PI = 3.141592653589793238462643383279502884;
#endif

typedef std::vector<double> dvec;

/**
 * Antenna patterns for different pulsars, which define how sensitive is Earth to the
 * pulsar oscillations
 *
 * @param extrinsic Intrinsic parameters of the GW source
 * @param pulsarUnitVector The unit vector into the direction of the pulsar
 */
dvec antennaPattern (dvec& intrinsic, dvec& pulsarUnitVector);

/**
 * Frequency evolution of the GW source
 *
 * @param t Time when the residual is evaluated
 * @param omega0 The initial frequency (at t=0)
 * @param M The chirp mass of the GW source
 */
inline double omegaReduced(const double omega0);

/**
 * Phase evolution of the GW source raised to the power of -1/3. This is just a
 * convenience, so that there are less computations done in some steps
 *
 * @param t Time when the residual is evaluated
 * @param omega0 The initial frequency (at t=0)
 * @param M The chirp mass of the GW source
 */
inline double Phi(const double t, const double omega0);

/**
 * Gravitational Wave Contribution to the residual:
 *
 * @param t Time when the residual is evaluated
 * @param intrinsic Intrinsic parameters of the GW source
 * @param omega The frequency of the source
 */
std::vector<double> gravWaveContrib (const double t, dvec& extrinsic, const double omega);

/**
 * Amplitude of the wave
 *
 * @param intrinsic Intrinsic parameters of the GW source
 */
dvec amplitude (dvec& intrinsic);

/**
 * The basis functions for the wave
 *
 * @param t Time when the residual is evaluated
 * @param extrinsic The extrinsic parameters of the source
 * @param pulsarUnitVector The unit vector into the direction of the pulsar
 */
dvec basis (const double t, dvec& intrinsic, dvec& pulsarUnitVector);

/**
 * The pulsar term for the residual
 *
 * @param t Time when the residual is evaluated
 * @param intrinsic Intrinsic parameters of the GW source
 * @param extrinsic The extrinsic parameters of the source
 * @param pulsarCoords The polar coordinates of a pulsar
 */
double pulsarTerm (const double t, dvec& extrinsic, dvec& intrinsic, dvec& pulsarCoords);

/**
 * Noise function
 *
 * @param t Time when the residual is evaluated
 * @param whiteNoise The amplitude of white noise in a pulsar
 * @param redNoise An array containing relevant red noise parameters
 * @param powerLawNoise An array containing relevant power law noise parameters
 */
double noise (const double t, const double whiteNoise, dvec& redNoise, dvec& powerLawNoise);

/**
 * The signal for one source excluding noise
 *
 * @param t Time when the residual is evaluated
 * @param params extrinsic and intrinsic vars in one array
 * @param L The distance to some selected pulsar
 * @param pulsarUnitVector The unit vector into the direction of the pulsar
 */
double individualSource (const double t, dvec& params, const double L, dvec& pulsarUnitVector);

/**
 * The residual for all the sources
 *
 * @param t Time when the residual is evaluated
 * @param N The number of sources
 * @param sources The parameters of the source(s). If there are more sources, then it is
 *        assumed, that they are put in rows of a matrix (a vector of a vector).
 * @param pulsarProperties The coordinates and the noise properties of the pulsar
 * @param noise whether to include noise in the signal
 */
double residual (const double t, const unsigned int N, std::vector<dvec>& sources, Pulsar& pulsarProperties, bool noise);
