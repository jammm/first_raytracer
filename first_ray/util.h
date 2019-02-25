#ifndef UTIL_H
#define UTIL_H

#ifdef _WIN32
#include <cstdlib>

float drand48()
{
    return float(rand()) / float(RAND_MAX);
}
#endif

#endif

