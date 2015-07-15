// -------------------------------------------------------------------------------
/*
    m4dUtilities.cpp

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dUtilities.h"
#include "m4dGlobalDefs.h"

namespace m4d {

// ---------------------------------------------------
/*!  Get system clock in micro-seconds.
 */
int64_t  get_system_clock() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (int64_t)tv.tv_sec * 1000000 + tv.tv_usec;
}


// ---------------------------------------------------
/*!  Tokenize file
 *  \param filename : name of file.
 *  \param tokens   : reference to vector of vector of string.
 *  \param useStandardIgnoreTokens : use standard ignore tokens ("#").
 *  \return  true : success.
 */
bool tokenizeFile(const std::string filename, std::vector<std::vector<std::string> > &tokens, bool useStandardIgnoreTokens) {
    std::ifstream in(filename.c_str());

    if (!in) {
        fprintf(stderr, "Cannot open file \"%s\"\n", filename.c_str());
        return false;
    }

    do {
        std::string line;
        getline(in, line);

        bool ignoreLine = false;
        if (!line.compare(0, 1, "")) {
            ignoreLine = true;
        }

        if (useStandardIgnoreTokens) {
            if (!line.compare(0, 1, "#")) {
                ignoreLine = true;
            }
        }

        if (!ignoreLine) {
            std::string buf;
            std::stringstream ss(line);

            std::vector<std::string> line_tokens;

            while (ss >> buf) {
                line_tokens.push_back(buf);
            }

            tokens.push_back(line_tokens);
        }
    } while (!in.eof());

    in.close();
    return true;
}

// ---------------------------------------------------
/*!  Tokenize file
 *  \param filename : name of file.
 *  \param ignores  : vector of strings which indicates lines that have to be ignored.
 *  \param tokens   : reference to vector of vector of string.
 *  \return  true : success.
 */
bool tokenizeFile(const std::string filename, const std::vector<std::string> ignores, std::vector<std::vector<std::string> > &tokens) {
    std::ifstream in(filename.c_str());

    if (!in) {
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        return false;
    }

    do {
        std::string line;
        getline(in, line);

        bool ignoreLine = false;

        for (unsigned int i = 0; i < ignores.size(); i++)
            if (!line.compare(0, ignores[i].length(), ignores[i])) {
                ignoreLine |= true;
            }

        if (!line.compare(0, 1, "")) {
            ignoreLine = true;
        }

        if (!ignoreLine) {
            std::vector<std::string> line_tokens;
            std::string buf;
            std::stringstream ss(line);

            while (ss >> buf) {
                line_tokens.push_back(buf);
            }

            tokens.push_back(line_tokens);
        }
    } while (!in.eof());

    in.close();
    return true;
}

// ---------------------------------------------------
/*!  get integer from tokens
 */
bool  getIntFromTokens(const std::vector<std::string> &tokenRow, std::string name, int &val) {
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = atoi(tokenRow[1].c_str());
        return true;
    }
    return false;
}

// ---------------------------------------------------
/*!  get integer array from tokens
 */
bool getIntFromTokensV(const std::vector<std::string> &tokenRow, std::string name, int num, int* val) {
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && (int)tokenRow.size() > (num + 1)) {
        for (int i = 0; i < num; i++) {
            val[i] = atoi(tokenRow[i + 1].c_str());
        }
        return true;
    }
    return false;
}

// ---------------------------------------------------
/*!  get double from tokens
 */
bool getDoubleFromTokens(const std::vector<std::string> &tokenRow, std::string name, double &val) {
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = atof(tokenRow[1].c_str());
        return true;
    }
    return false;
}

// ---------------------------------------------------
/*!  get double array from tokens
 */
bool getDoubleFromTokensV(const std::vector<std::string> &tokenRow, std::string name, int num, double* val) {
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && (int)tokenRow.size() > (num + 1)) {
        for (int i = 0; i < num; i++) {
            val[i] = atof(tokenRow[i + 1].c_str());
        }
        return true;
    }
    return false;
}

// ---------------------------------------------------
/*!  get string from tokens
 */
bool getStringFromTokens(const std::vector<std::string> &tokenRow, std::string name, std::string &val) {
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = tokenRow[1];
        return true;
    }
    return false;
}

/*! Write a binary float array.
 *   \param   filename  :  name of the file.
 *   \param      array  :  pointer to float array.
 *   \param          x  :  width of array.
 *   \param          y  :  height of array.
 *   \param          c  :  number of channels.
 *   \return      true  :  success.
 *   \return     false  :  error occured.
 */
bool  writeFloatArray(std::string filename, const float* array, int x, int y, int c) {
    std::ofstream out(filename.c_str(), std::ios::binary);
    if (!out.is_open()) {
        fprintf(stderr, "Can't open file %s for output.\n", filename.c_str());
        return false;
    }

    std::string hdr = "HEAD";
    char buf[5];
#ifdef _WIN32
    sprintf_s(buf, "%s", hdr.c_str());
#else
    sprintf(buf, "%s", hdr.c_str());
#endif
    std::string hdrd = "DATA";
    char bufd[5];
#ifdef _WIN32
    sprintf_s(bufd, "%s", hdrd.c_str());
#else
    sprintf(bufd, "%s", hdrd.c_str());
#endif

    int hdrSize = 24;
    out.write((char*)&buf[0], sizeof(char) * 4);
    out.write((char*)&hdrSize, sizeof(int));
    out.write((char*)&x, sizeof(int));
    out.write((char*)&y, sizeof(int));
    out.write((char*)&c, sizeof(int));
    out.write((char*)&bufd[0], sizeof(char) * 4);
    out.write((char*)&array[0], sizeof(float)*x * y * c);
    out.close();
    return true;
}

/*! Read a binary float array.
 *   \param   filename  :  name of the file.
 *   \param          x  :  reference to width of array.
 *   \param          y  :  reference to height of array.
 *   \param          c  :  reference to number of channels.
 *   \return      true  :  success.
 *   \return     false  :  error occured.
 */
float* readFloatArray(std::string filename, int &x, int &y, int &c) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.is_open()) {
        fprintf(stderr, "Can't open file %s for reading.\n", filename.c_str());
        return NULL;
    }

    char buf[5];
    char bufd[5];
    for (int i = 0; i < 4; i++) {
        buf[i] = in.peek();
        in.seekg(i + 1);
    }

    float* array = NULL;
    if (strncmp(buf, "HEAD", 4) == 0) {
        int hdrSize;
        in.read((char*)&hdrSize, sizeof(int));  // header size
        in.read((char*)&x, sizeof(int));        // resX
        in.read((char*)&y, sizeof(int));        // resY
        in.read((char*)&c, sizeof(int));        // numChannels
        in.read((char*)&bufd, sizeof(char) * 4);
        if (strncmp(bufd, "DATA", 4) == 0) {
            if (array != NULL) {
                delete [] array;
            }
            array = new float[x * y * c];
            std::streamsize size = x * y * c * sizeof(float);
            in.read((char*)&array[0], size);
        }
        in.close();
    }

    return array;
}

double radians(double phi) {
    return phi * DEG_TO_RAD;
}

double degree(double phi) {
    return phi * RAD_TO_DEG;
}


int find_nat_tetrad_type( const char* name) {
    unsigned int n = 0;
    while(n < NUM_ENUM_NAT_TETRAD_TYPES) {
        if (strcmp(name,stl_nat_tetrad_types[n])==0) {
            return n;
        }
        n++;
    }
    return -1;
}

} // end namespace m4d

