/**
 * @file    m4dPlatform.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
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
