/**
 * @file    m4dUtilities.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dUtilities.h"
#include "m4dGlobalDefs.h"
#ifdef _WIN32
#include <windows.h>
#endif

namespace m4d {

int64_t get_system_clock()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return static_cast<int64_t>(tv.tv_sec) * 1000000 + tv.tv_usec;
}

bool tokenizeFile(
    const std::string& filename, std::vector<std::vector<std::string>>& tokens, bool useStandardIgnoreTokens)
{
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
bool tokenizeFile(
    const std::string& filename, const std::vector<std::string>& ignores, std::vector<std::vector<std::string>>& tokens)
{
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
bool getIntFromTokens(const std::vector<std::string>& tokenRow, std::string name, int& val)
{
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = atoi(tokenRow[1].c_str());
        return true;
    }
    return false;
}

bool getIntFromTokensV(const std::vector<std::string>& tokenRow, std::string name, int num, int* val)
{
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && (int)tokenRow.size() > (num + 1)) {
        for (int i = 0; i < num; i++) {
            val[i] = atoi(tokenRow[i + 1].c_str());
        }
        return true;
    }
    return false;
}

bool getDoubleFromTokens(const std::vector<std::string>& tokenRow, std::string name, double& val)
{
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = atof(tokenRow[1].c_str());
        return true;
    }
    return false;
}

bool getDoubleFromTokensV(const std::vector<std::string>& tokenRow, std::string name, int num, double* val)
{
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && static_cast<int>(tokenRow.size()) > (num + 1)) {
        for (int i = 0; i < num; i++) {
            val[i] = atof(tokenRow[i + 1].c_str());
        }
        return true;
    }
    return false;
}

bool getStringFromTokens(const std::vector<std::string>& tokenRow, std::string name, std::string& val)
{
    std::string baseString = tokenRow[0];
    if (baseString.compare(name) == 0 && tokenRow.size() > 1) {
        val = tokenRow[1];
        return true;
    }
    return false;
}

bool writeFloatArray(std::string filename, const float* array, int x, int y, int c)
{
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
    out.write((char*)&array[0], sizeof(float) * x * y * c);
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
float* readFloatArray(std::string filename, int& x, int& y, int& c)
{
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.is_open()) {
        fprintf(stderr, "Can't open file %s for reading.\n", filename.c_str());
        return nullptr;
    }

    char buf[5];
    char bufd[5];
    for (int i = 0; i < 4; i++) {
        buf[i] = in.peek();
        in.seekg(i + 1);
    }

    float* array = nullptr;
    if (strncmp(buf, "HEAD", 4) == 0) {
        int hdrSize;
        in.read((char*)&hdrSize, sizeof(int)); // header size
        in.read((char*)&x, sizeof(int)); // resX
        in.read((char*)&y, sizeof(int)); // resY
        in.read((char*)&c, sizeof(int)); // numChannels
        in.read((char*)&bufd, sizeof(char) * 4);
        if (strncmp(bufd, "DATA", 4) == 0) {
            if (array != nullptr) {
                delete[] array;
            }
            array = new float[x * y * c];
            std::streamsize size = x * y * c * sizeof(float);
            in.read((char*)&array[0], size);
        }
        in.close();
    }

    return array;
}

double M4D_CALL radians(double phi)
{
    return phi * DEG_TO_RAD;
}

double M4D_CALL degree(double phi)
{
    return phi * RAD_TO_DEG;
}

int M4D_CALL find_nat_tetrad_type(const char* name)
{
    unsigned int n = 0;
    while (n < NUM_ENUM_NAT_TETRAD_TYPES) {
        if (strcmp(name, stl_nat_tetrad_types[n]) == 0) {
            return n;
        }
        n++;
    }
    return -1;
}

bool CopyString(const char* src, char*& dest)
{
    if (src == nullptr) {
        return false;
    }

    if (dest != nullptr) {
        delete[] dest;
    }
    size_t len = strlen(src);

    bool isOkay = true;
#ifdef _WIN32
    dest = new char[len + 4];
    isOkay &= (strncpy_s(dest, len + 4, src, len) == 0);
#else
    dest = new char[len + 2];
    isOkay &= (strcpy(dest, src) != nullptr);
#endif
    return isOkay;
}

} // end namespace m4d
