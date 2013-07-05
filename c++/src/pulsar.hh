#include <vector>

typedef std::vector<double> dvec;

#ifndef _PULSAR_HXX_
#define _PULSAR_HXX_

class Pulsar {
    private:
        double mDistance,
               mTheta,
               mPhi,
               mWhiteNoise;
        dvec mRedNoise,
             mPowerLawNoise,
             mTimes;

    public:
        Pulsar ();

        ~Pulsar () {}

        // Setters
        
        /**
         * Set the angular coordinates of a pulsar.
         *
         * @param theta Azimuthal angle \theta
         * @param phi Polar angle \phi
         */
        void setAngles (double theta, double phi) { mTheta = theta, mPhi = phi; }

        /**
         * Set the distance of the pulsar
         *
         * This has to be separated from the other coordinates, as we might want to vary
         * when we are mathching the pulsar term
         *
         * @param x The distance from Solar Barycentre to the Pulsar
         */
        void setDistance (double x) { mDistance = x; }

        /**
         * Set the white noise amplitude
         *
         * @param N amplitude of the fluctuations
         */
        void setWhiteNoise (double N) { mWhiteNoise = N; }

        /**
         * Set the red noise parameters
         *
         * @param N amplitude of the fluctuations
         * @param f specific frequency of the noise
         */
        void setRedNoise (double N, double f);

        /**
         * Set the power law noise parameters
         *
         * @param N amplitude of the fluctuations
         * @param gamma spectral index
         */
        void setPowerLawNoise (double N, double gamma);

        // Getters
        /**
         * Get a White noise amplitude
         */
        double getWhiteNoise () { return mWhiteNoise; }

        /**
         * Get red noise parameters
         */
        dvec getRedNoise () { return mRedNoise; }

        /**
         * Get power law noise parameters
         */
        dvec getPowerLawNoise () { return mPowerLawNoise; }

        /**
         * Get a distance for the pulsar
         */
        double getDistance () { return mDistance; }

        /**
         * Get a unit vector for the pulsar
         */
        dvec getUnitVector ();
};

#endif

namespace pulsarGrid {
    /**
     * Randomize the pulsar grid
     *
     * @param pulsars Is the vector template for generating the data in
     * @param N The number of pulsar to be put
     * @param range The ranges of coordinates. It should be an array with 6 entries
     * @param whiteNoise White Noise amplitude for a pulsar
     */
    void randomizeData (std::vector<Pulsar>& Grid, unsigned int N, dvec range, double wnoise);

    /**
     * Sampling schedule generation
     *
     * @param pulsars Is the vector template for generating the data in
     * @param initialTimes A list of times for a first measurement to be taken. The array
     *  index is the same as a pulsar index.
     * @param tFinal Finishing time
     * @param tMin A minimum time break between 2 samples.
     * @param tMax A maximum time break between 2 samples.
     * @param indices_out The output array for the indices.
     * @param Times_out The output array for the time schedule.
     */
    void generateSchedule (std::vector<double>& initialTimes, 
                           double tFinal, double tMin, double tMax,
                           std::vector<unsigned short>& indices_out, 
                           std::vector<double>& Times_out);
}
