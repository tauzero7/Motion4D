/**
 * @file    m4dPlatform.h
 * @author  Thomas Mueller
 *
 * @brief  Platform dependend functions
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_PLATFORM_H
#define M4D_PLATFORM_H

#ifdef _WIN32
#include <time.h>
//#include <windows.h>
#include <stdio.h>
typedef signed long long int int64_t;
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL
#endif

struct timezone {
    int tz_minuteswest; /* minutes W of Greenwich */
    int tz_dsttime; /* type of dst correction */
};

/*! Windows version of Linux gettimeofday
 */
int gettimeofday(struct timeval* tv, struct timezone* tz);

#else
#include <sys/time.h>
#endif //_WIN32

#endif
