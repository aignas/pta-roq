#include <vector>

typedef std::vector<double> dvec

class PulsarGrid {
    private:
        std::vector<dvec> mAngle,
                          mRedNoise,
                          mPowerLawNoise;
        dvec mDistance;
        dvec mWhiteNoise;
    public:
        // Constructors and destructors
        PulsarGrid (unsigned int, dvec, double);

        PulsarGrid ();

        ~PulsarGrid () {}

        // Setters
        //
        // This randomizes the data, the default option
        void randomizeData (unsigned int, dvec, dvec);

        void setSize (unsigned int);

        void setWhiteNoise (dvec);

        void setAngles (std::vector<dvec>);

        // Getters
        /**
         * Get a White noise amplitude
         *
         * @param idx The label of pulsar
         */
        double getWhiteNoise (unsigned int idx) { return mWhiteNoise[idx]; }

        /**
         * Get red noise parameters
         *
         * @param idx The label of pulsar
         */
        dvec getRedNoise (unsigned int idx) { return mRedNoise[idx]; }

        /**
         * Get power law noise parameters
         *
         * @param idx The label of pulsar
         */
        dvec getPowerLawNoise (unsigned int idx) { return mPowerLawNoise[idx]; }

        /**
         * Get a distance for the pulsar
         *
         * @param idx The label of pulsar
         */
        double getDistance (unsigned int idx) { return mDistance[idx]; }

        /**
         * Get the number of pulsars
         */
        unsigned int getNumber () { return mAngle.size(); }

        /**
         * Get a unit vector for the pulsar
         *
         * @param idx The label of pulsar
         */
        dvec getUnitVector (unsigned int);
};
