#include <string>
#include <vector>

#include "pulsar.hh"

/**
 * Convert the csv file into an array of doubles
 */
int csv2array (std::string filename, std::vector<double> & array);

/**
 * Read pulsar data
 *
 * @param filename A filename to read the pulsar data from
 * @param pulsars The array to store pulsar array in
 * @param The delimiter in the file
 */
int csv2pulsar (const std::string filename, std::vector<Pulsar> & pulsars, const std::string delim);

/**
 * write pulsar data
 *
 * @param filename A filename to write the pulsar data to
 * @param pulsars The array to use for writing
 * @param The delimiter in the file
 */
int pulsar2csv (const std::string filename, std::vector<Pulsar> & pulsars, const std::string delim);

/**
 * write the schedule data
 *
 * @param filename A filename to write the pulsar data to
 * @param iarray An integer array
 * @param darray The double array
 * @param The delimiter in the file
 */
int arraysShortDouble2csv (const std::string filename, std::vector<unsigned short> & iarray, std::vector<double> & darray, const std::string delim);
int arraysLongDouble2csv (const std::string filename, std::vector<long> & iarray, std::vector<double> & darray, const std::string delim);

/**
 * read the schedule data
 *
 * @param filename A filename to write the pulsar data to
 * @param iarray The integer array
 * @param darray The double array
 * @param The delimiter in the file
 */
int csv2arraysShortDouble (const std::string filename, std::vector<unsigned short> & iarray, std::vector<double> & darray, const std::string delim);
int csv2arraysLongDouble (const std::string filename, std::vector<long> & iarray, std::vector<double> & darray, const std::string delim);

int csv2paramRanges (const std::string filename, std::vector<double> & params_min,
        std::vector<double> & params_max, std::vector<unsigned int> & params_N, const std::string delim);

/**
 * write the source data
 *
 * @param filename A filename to write the pulsar data to
 * @param sources An array of sources
 * @param The delimiter in the file
 */
int arrayArrayDouble2csv (const std::string filename, std::vector<std::vector<double> > & sources, const std::string delim);


/**
 * read the source data
 *
 * @param filename A filename to write the pulsar data to
 * @param sources An array of sources
 * @param The delimiter in the file
 */
int csv2arrayArrayDouble (const std::string filename, std::vector<std::vector<double> > & sources, const std::string delim);
int csv2arrayArrayUnsignedInt (const std::string filename, std::vector<std::vector<unsigned int> > & sources, const std::string delim);

int arrayDouble2csv (const std::string filename, std::vector<double> & array, const std::string delim);

int csv2arrayDouble (const std::string filename, std::vector<double> & array, const std::string delim);

int parseDataRC (const std::string & filename, std::string & stamp, std::vector<std::string> & fnames, std::string & delim);

int parseRoqRC (const std::string & filename, std::string & stamp, std::string & in_stamp, std::vector<std::string> & fnames, std::string & delim);

/**
 * Trim a string from the left
 *
 * @param s string to trim from the left
 */
namespace helper {
    std::string &ltrim(std::string &s);

    /**
     * Trim a string from the right
     *
     * @param s string to trim
     */
    std::string &rtrim(std::string &s);

    /**
     * Trim a string from the both sides
     *
     * @param s string to trim
     */
    std::string &trim(std::string &s);

    /**
     * Split a string into a vector of strings at the place of some delimiter
     *
     * @param s string to split
     * @param delim the delimiter
     * @param elems the output vector
     */
    void split(const std::string &s, std::string delim, std::vector<std::string> &elems);

    /**
     * Convert a string to double
     *
     * @param s string to convert
     */
    double convertToDouble(std::string const& s);

    /**
     * Convert a string to unsigned
     *
     * @param s string to convert
     */
    unsigned int convertToUnsignedInt (std::string const& s);

    /**
     * Convert a string to unsigned long
     *
     * @param s string to convert
     */
    unsigned long convertToUnsignedLong (std::string const& s);

    /**
     * Convert a string to unsigned short
     *
     * @param s string to convert
     */
    unsigned short convertToUnsignedShort (std::string const& s);

    int removeComments (std::ifstream & fin, std::string & value);

    std::string getKey (std::ifstream & fin, const char * key);
}
