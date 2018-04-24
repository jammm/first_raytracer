#ifndef UTIL_H
#define UTIL_H

#include <cstdlib>

float drand48()
{
    return float(rand()) / float(RAND_MAX);
}

#endif

