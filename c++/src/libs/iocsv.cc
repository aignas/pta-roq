#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

#include "pulsar.hh"

namespace helper {
    // trim from start
    std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                    std::not1(std::ptr_fun<int, int>(std::isspace))));

        return s;
    }

    // trim from end
    std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                    std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
                s.end());

        return s;
    }

    // trim from both ends
    std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
    }

    void split(const std::string &s, std::string delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;

        // Clear the vector
        elems.clear();

        while (std::getline(ss, item, delim[0])) {
            elems.push_back(item);
        }
    }

    double convertToDouble(std::string const& s) {
        std::istringstream i(s);
        double x;
        if (!(i >> x)) {
            return 0;
        }

        return x;
    }

    unsigned int convertToUnsignedInt (std::string const& s) {
        std::istringstream i(s);
        unsigned int x;
        if (!(i >> x)) {
            return 0;
        }

        return x;
    }

    unsigned long convertToUnsignedLong (std::string const& s) {
        std::istringstream i(s);
        unsigned long x;
        if (!(i >> x)) {
            return 0;
        }

        return x;
    }

    unsigned short convertToUnsignedShort (std::string const& s) {
        std::istringstream i(s);
        unsigned short x;
        if (!(i >> x)) {
            return 0;
        }

        return x;
    }

    int removeComments (std::ifstream & fin, std::string & value) {
        bool GotIt = false;

        if (fin.is_open()) {
            while (not GotIt) {
                std::getline(fin, value);

                if (value[0] != '#' and not (value[0] == '/' and value[1] == '/')) {
                    GotIt = true;
                }
            }
        } else {
            std::cerr << "Error opening file." << std::endl;
            return 1;
        }

        return 0;
    }

    std::string getKey (std::ifstream & fin, const char * key) {
        std::string comp = key;
        bool GotIt = false;
        std::string value;
        std::vector<std::string> pair;

        if (fin.is_open()) {
            while (not GotIt and fin.good()) {
                std::getline(fin, value);

                helper::removeComments(fin, value);
                helper::split(value, ":", pair);
                helper::trim(pair[0]);
                helper::trim(pair[1]);

                if (comp.compare(pair[0])) {
                    GotIt = true;
                }
            }

            fin.close();
        }

        if (GotIt) {
            return pair[1];
        } else {
            std::cerr << "A key \"" << key
                << "\" was not found, please check your configuration file"
                << std::endl;
            return std::string();
        }
    }
}

int csv2arrayDouble (const std::string filename, std::vector<double> array, const std::string delim) {
    std::ifstream fin (filename);
    std::string value;

    if (fin.is_open()) {
        std::cout << "Reading array data from " << filename << std::endl;

        array.clear();

        while (fin.good()) {
            std::getline(fin, value, delim[0]);

            array.push_back(helper::convertToDouble(value));
        }

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int arrayDouble2csv (const std::string filename, std::vector<double> array, const std::string delim) {
    std::ofstream fout (filename);

    if (fout.is_open()) {
        std::cout << "Saving array data into " << filename << std::endl;

        fout << array.at(0);

        for (unsigned i = 1; i < array.size(); i++) {
            fout << delim << array.at(i);
        }

        fout.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int pulsar2csv (const std::string filename, std::vector<Pulsar> & pulsars, const std::string delim) {
    std::ofstream fout (filename);
    const unsigned int pulsarNumber = pulsars.size();

    if (fout.is_open()) {
        std::cout << "Saving pulsar data into " << filename << std::endl;
        std::vector<double> tmp;
        for (unsigned i = 0; i < pulsarNumber; i++) {
            pulsars[i].getAll(tmp);

            // Have delimters separating the entries
            for (unsigned j = 0; j < tmp.size() - 1; j++) {
                fout << tmp.at(j) << delim;
            }

            fout << tmp.back();

            if (i != pulsarNumber -1) {
                fout << std::endl;
            }
        }
        fout.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int csv2pulsar (const std::string filename, std::vector<Pulsar> & pulsars, const std::string delim) {
    std::ifstream fin (filename);
    std::string value;

    // check if the file is really open
    if (fin.is_open()) {
        std::cout << "Reading pulsar data from " << filename << std::endl;

        std::vector<std::string> elements;
        std::vector<double> elements_dbl;
        pulsars.clear();

        while (fin.good()) {
            std::getline(fin, value);
            helper::split (value, delim, elements);

            // Cast everithing into doubles
            elements_dbl.resize(elements.size());
            for (unsigned i = 0; i < elements.size(); i++) {
                elements_dbl[i] = helper::convertToDouble(elements[i]);
            }

            // Add a pulsar to the list
            pulsars.push_back(Pulsar (elements_dbl));
        }

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int arraysLongDouble2csv (const std::string filename, std::vector<long> & iarray, std::vector<double> & darray, const std::string delim) {
    std::ofstream fout (filename);

    if (fout.is_open()) {
        std::cout << "Saving time-stamp data into " << filename << std::endl;
        for (unsigned i = 0; i < iarray.size(); i++) {
            fout << iarray.at(i) << delim << darray.at(i);

            if (i != iarray.size() - 1) {
                fout << std::endl;
            }
        }
        fout.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int arraysShortDouble2csv (const std::string filename, std::vector<unsigned short> & iarray, std::vector<double> & darray, const std::string delim) {
    std::ofstream fout (filename);

    if (fout.is_open()) {
        std::cout << "Saving time-stamp data into " << filename << std::endl;
        for (unsigned i = 0; i < iarray.size(); i++) {
            fout << iarray.at(i) << delim << darray.at(i);

            if (i != iarray.size() - 1) {
                fout << std::endl;
            }
        }
        fout.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int arrayArrayDouble2csv (const std::string filename, std::vector<std::vector<double> > & array, const std::string delim) {
    std::ofstream fout (filename);

    if (fout.is_open()) {
        std::cout << "Saving array data into " << filename << std::endl;
        for (unsigned i = 0; i < array.size(); i++) {
            fout << array.at(i).at(0);

            for (unsigned j = 1; j < array.at(i).size(); j++) {
                fout << delim << array.at(i).at(j);
            }

            if (i != array.size() - 1) {
                fout << std::endl;
            }
        }

        fout.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int csv2arraysShortDouble (const std::string filename, std::vector<unsigned short> & iarray, std::vector<double> & darray, const std::string delim) {
    std::ifstream fin (filename);
    std::string value;

    // check if the file is really open
    if (fin.is_open()) {
        std::cout << "Reading time-stamp data from " << filename << std::endl;

        iarray.clear();
        darray.clear();

        while (fin.good()) {
            // Index
            std::getline(fin, value, delim[0]);
            iarray.push_back(helper::convertToUnsignedShort(value));

            std::getline(fin, value, delim[0]);
            darray.push_back(helper::convertToDouble(value));
        }

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int csv2arrayArrayDouble (const std::string filename, std::vector<std::vector<double> > & array, const std::string delim) {
    std::ifstream fin (filename);
    std::string value;

    // check if the file is really open
    if (fin.is_open()) {
        std::cout << "Reading array data from " << filename << std::endl;

        std::vector<std::string> elements;
        std::vector<double> elements_dbl;
        array.clear();

        while (fin.good()) {
            std::getline(fin, value);
            helper::split (value, delim, elements);

            // Cast everithing into doubles
            elements_dbl.resize(elements.size());
            for (unsigned i = 0; i < elements.size(); i++) {
                elements_dbl[i] = helper::convertToDouble(elements[i]);
            }

            // Add a pulsar to the list
            array.push_back(elements_dbl);
        }

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    return 0;
}

int parseRoqRC (const std::string & filename, std::string & stamp, std::string & in_stamp, std::vector<std::string> & fnames, std::string & delim) {
    std::string ext, value, separator;
    std::vector<std::string> dir, prefix;
    std::ifstream fin (filename);

    // check if the file is really open
    if (fin.is_open()) {
        std::cout << "Reading program parameters from " << filename << std::endl;

        dir.push_back(helper::getKey (fin, "dir.out"));
        dir.push_back(helper::getKey (fin, "dir.in"));
        prefix.push_back(helper::getKey (fin, "prefix.out"));
        prefix.push_back(helper::getKey (fin, "prefix.in"));

        ext = helper::getKey (fin, "ext");
        separator = helper::getKey (fin, "sep");
        delim = helper::getKey (fin,"delim");

        // This ensures no bugs in the future
        fnames.clear();
        fnames.push_back(helper::getKey(fin, "infile.pulsar"));
        fnames.push_back(helper::getKey(fin, "infile.sched"));
        fnames.push_back(helper::getKey(fin, "infile.resid"));

        fnames.push_back(helper::getKey(fin, "outfile.sched"));
        fnames.push_back(helper::getKey(fin, "outfile.eim"));
        fnames.push_back(helper::getKey(fin, "outfile.resid"));
        fnames.push_back(helper::getKey(fin, "outfile.rb"));
        fnames.push_back(helper::getKey(fin, "outfile.rbparams"));

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    std::stringstream fnameFull;
    // Generate the input filenames
    for (unsigned i = 0; i < 3; i++) {
        fnameFull.str(std::string());
        fnameFull << dir[1] << "/" << prefix[1] << separator << in_stamp << separator << fnames[i] << ext;
        fnames[i] = fnameFull.str();
    }

    // Generate the filenames
    for (unsigned i = 3; i < fnames.size(); i++) {
        fnameFull.str(std::string());
        fnameFull << dir[0] << "/" << prefix[0] << separator << stamp << separator << fnames[i] << ext;
        fnames[i] = fnameFull.str();
    }

    return 0;
}

int parseDataRC (const std::string & filename, std::string & stamp, std::vector<std::string> & fnames, std::string & delim) {
    std::string ext, value, separator;
    std::vector<std::string> dir, prefix;
    std::ifstream fin (filename);

    // check if the file is really open
    if (fin.is_open()) {
        std::cout << "Reading program parameters from " << filename << std::endl;

        dir.push_back(helper::getKey (fin, "dir.out"));
        dir.push_back(helper::getKey (fin, "dir.in"));
        prefix.push_back(helper::getKey (fin, "prefix.out"));
        prefix.push_back(helper::getKey (fin, "prefix.in"));

        ext = helper::getKey (fin, "ext");
        separator = helper::getKey (fin, "sep");
        delim = helper::getKey (fin,"delim");

        // This ensures no bugs in the future
        fnames.clear();
        fnames.push_back(helper::getKey(fin, "infile.source"));
        fnames.push_back(helper::getKey(fin, "outfile.pulsar"));
        fnames.push_back(helper::getKey(fin, "outfile.sched"));
        fnames.push_back(helper::getKey(fin, "outfile.resid"));

        fin.close();
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1;
    }

    std::stringstream fnameFull;
    // Generate the filenames for output data
    for (unsigned i = 0; i < 1; i++) {
        fnameFull.str(std::string());
        fnameFull << dir[1] << "/" << prefix[1] << separator << fnames[i] << ext;
        fnames[i] = fnameFull.str();
    }

    // Generate the filenames for output data
    for (unsigned i = 1; i < fnames.size(); i++) {
        fnameFull.str(std::string());
        fnameFull << dir[0] << "/" << prefix[0] << separator << stamp << separator << fnames[i] << ext;
        fnames[i] = fnameFull.str();
    }

    return 0;
}