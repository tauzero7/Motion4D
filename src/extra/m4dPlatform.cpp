// -------------------------------------------------------------------------------
/*
    m4dPlatform.cpp

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

#include "m4dPlatform.h"

#ifdef _WIN32
#include <windows.h>

// ---------------------------------------------------
//    gettimeofday
// ---------------------------------------------------
int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (tv != nullptr) {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres /= 10; /*convert into microseconds*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (tz != nullptr) {
        if (!tzflag) {
            _tzset();
            tzflag++;
        }
#ifndef _WIN32
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
#else
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
        // tz->tz_minuteswest = _get_timezone / 60;
        // tz->tz_dsttime = _get_daylight;
#endif
    }

    return 0;
}
#endif
