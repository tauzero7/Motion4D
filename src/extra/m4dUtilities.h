// --------------------------------------------------------------------------------
/*
    m4dUtilities.h

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

/*!  \file   m4dUtilities.h
     \brief  Utility functions.

*/
// --------------------------------------------------------------------------------

#ifndef M4D_UTILITIES_H
#define M4D_UTILITIES_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cstring>
#include <sys/types.h>

#include "m4dPlatform.h"

namespace m4d {


int64_t  get_system_clock();


bool tokenizeFile(const std::string filename, std::vector<std::vector<std::string> > &tokens,
                  bool useStandardIgnoreTokens = true);

bool tokenizeFile(const std::string filename, const std::vector<std::string> ignores,
                  std::vector<std::vector<std::string> > &tokens);

bool getIntFromTokens(const std::vector<std::string> &tokenRow, std::string name, int &val);
bool getIntFromTokensV(const std::vector<std::string> &tokenRow, std::string name, int num, int* val);
bool getDoubleFromTokens(const std::vector<std::string> &tokenRow, std::string name, double &val);
bool getDoubleFromTokensV(const std::vector<std::string> &tokenRow, std::string name, int num, double* val);
bool getStringFromTokens(const std::vector<std::string> &tokenRow, std::string name, std::string &val);

bool     writeFloatArray(std::string filename, const float* array, int x, int y, int c);
float*   readFloatArray(std::string filename, int &x, int &y, int &c);

} // end namespace m4d

#endif


























