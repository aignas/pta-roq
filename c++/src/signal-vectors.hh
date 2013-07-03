/* This is a major rewrite of the code to make it more parallel
 *
 */

#include <vector>

class UnitVectors {
    // We have two angles and evaluate vectors on the fly
    private:
        double mTheta, mPhi;

    public:
        UnitVectors (double, double);

        // Empty constructor
        UnitVectors ();

        // Desctructor
        ~UnitVectors () {}

        std::vector<double> Omega ();
        std::vector<double> m ();
        std::vector<double> n ();
};

