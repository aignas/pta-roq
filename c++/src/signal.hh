// Header file for signal generation

class PulsarGrid {
    private:
        std::vector<cav::Vec3> mR;
        std::vector<std::vector<double>> mNoise;
    public:
        PulsarGrid (int N, std::vector<double> ranges, std::vector<double> noise);

        std::vector<double>  (int idx);
        double getWhiteNoise (int idx);
        double getRedNoise   (int idx);
        double getPowLawNoise(int idx);
        std::vector<double> getAngles (int idx);
        cav::Vec3 getUnitVector (int idx);
        
        // Get the distance of the pulsars
        double getLength (int idx) { return mR[idx].x(); }

        // Get the number of pulsars
        int getNumber () { return mR.size(); }
}
