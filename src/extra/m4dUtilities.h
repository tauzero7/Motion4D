/**
 * @file    m4dUtilities.h
 * @author  Thomas Mueller
 *
 * @brief  Utility functions.
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_UTILITIES_H
#define M4D_UTILITIES_H

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <vector>

#include "m4dPlatform.h"

#if defined _WIN32 && defined(M4D_LIB)
#if defined(m4d_EXPORTS) || defined(m4d_lua_EXPORTS) || defined(m4dd_EXPORTS) || defined(m4d_luad_EXPORTS)             \
    || defined(_m4d_EXPORTS)
#define API_M4D_EXPORT __declspec(dllexport)
#else
#define API_M4D_EXPORT __declspec(dllimport)
#endif
#define M4D_CALL __stdcall
#else // _WIN32
#define M4D_CALL
#define API_M4D_EXPORT
#endif // _WIN32

namespace m4d {

/// Get system clock in micro-seconds.
int64_t get_system_clock();

/*!  Tokenize file
 *  \param filename : name of file.
 *  \param tokens   : reference to vector of vector of string.
 *  \param useStandardIgnoreTokens : use standard ignore tokens ("#").
 *  \return  true : success.
 */
bool tokenizeFile(
    const std::string& filename, std::vector<std::vector<std::string>>& tokens, bool useStandardIgnoreTokens = true);

bool tokenizeFile(const std::string& filename, const std::vector<std::string>& ignores,
    std::vector<std::vector<std::string>>& tokens);

bool getIntFromTokens(const std::vector<std::string>& tokenRow, std::string name, int& val);
bool getIntFromTokensV(const std::vector<std::string>& tokenRow, std::string name, int num, int* val);
bool getDoubleFromTokens(const std::vector<std::string>& tokenRow, std::string name, double& val);
bool getDoubleFromTokensV(const std::vector<std::string>& tokenRow, std::string name, int num, double* val);
bool getStringFromTokens(const std::vector<std::string>& tokenRow, std::string name, std::string& val);

/*! Write a binary float array.
 *   \param   filename  :  name of the file.
 *   \param      array  :  pointer to float array.
 *   \param          x  :  width of array.
 *   \param          y  :  height of array.
 *   \param          c  :  number of channels.
 *   \return      true  :  success.
 *   \return     false  :  error occured.
 */
bool writeFloatArray(std::string filename, const float* array, int x, int y, int c);
float* readFloatArray(std::string filename, int& x, int& y, int& c);

// prototype of functions
double API_M4D_EXPORT radians(double phi);
API_M4D_EXPORT double degree(double phi);

int API_M4D_EXPORT find_nat_tetrad_type(const char* name);

bool CopyString(const char* src, char*& dest);

/**
 * @brief Safely delete 1D arrays generated with 'new'.
 *
 *  Each pointer to a 1D array should be immediately
 *  initialized or set to nullptr. This method checks
 *  if pointer is not 'nullptr', deletes the array,
 *  and sets the pointer to 'nullptr'.
 *
 * @tparam T  Pointer to data array.
 */
template <typename T> void SafeDelete(T*& ptr)
{
    if (ptr != nullptr) {
        delete[] ptr;
        ptr = nullptr;
    }
}

} // end namespace m4d

#endif
