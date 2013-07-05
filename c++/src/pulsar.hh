#include <vector>

typedef std::vector<double> dvec;

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

        /**
         * Set the schedule
         *
         * @param t An array containing all the time stamps
         */
        void setSchedule (dvec t) { mTimes = t; }

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

        /**
         * Get the schedule
         */
        dvec getSchedule () { return mTimes; }
};

namespace PulsarGrid {
        /**
         * Randomize the pulsar grid
         *
         * @param N The number of pulsar to be put
         * @param range The ranges of coordinates. It should be an array with 6 entries
         * @param whiteNoise White Noise amplitude for a pulsar
         */
        void randomizeData (unsigned int N, dvec range, double whiteNoise);

        /**
         * Sampling schedule generation
         *
         * @param initialTimes A list of times for a first measurement to be taken. The array
         *  index is the same as a pulsar index.
         * @param tFinal Finishing time
         * @param tMin A minimum time break between 2 samples.
         * @param tMax A maximum time break between 2 samples.
         */
        void generateSchedule (std::vector<double> initialTimes, double tFinal, double tMin, double tMax);
};
