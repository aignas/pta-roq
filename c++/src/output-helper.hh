/**
 * This will contain all the helpers for array and other output as it is useful for
 * debugging
 */

#include <string>
#include <vector>

/**
 * Output a vector to the stdout
 *
 * @param A the vector to output
 * @param separator A separator which will separate values
 */
void outputVector (std::vector<double> & A, std::string separator);
void outputVector (std::vector<int> & A, std::string separator);

/**
 * Output a square matrix to the stdout
 *
 * @param A the matrix to output
 * @param sep1 A separator which will separate columns
 * @param sep2 A separator which will separate rows
 */
void outputMatrix (std::vector<double> & A, std::string sep1, std::string sep2);
void outputMatrix (std::vector<int> & A, std::string sep1, std::string sep2);
