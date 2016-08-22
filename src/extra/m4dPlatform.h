// --------------------------------------------------------------------------------
/*
    m4dPlatform.h

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

/*!  \file   m4dPlatform.h
     \brief  Platform dependend functions

*/
// --------------------------------------------------------------------------------

#ifndef M4D_PLATFORM_H
#define M4D_PLATFORM_H

#ifdef _WIN32
#include <time.h>
//#include <windows.h>
#include <stdio.h>
typedef signed long long int int64_t;
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone {
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

/*! Windows version of Linux gettimeofday
 */
int gettimeofday(struct timeval *tv, struct timezone *tz);

#else
#include <sys/time.h>
#endif //_WIN32


#endif
