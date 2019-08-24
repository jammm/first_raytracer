#ifndef SH_SAMPLE_H
#define SH_SAMPLE_H

#include "geometry.h"
#include "util.h"

namespace PRT
{

    constexpr int n_bands = 7;
    constexpr int n_coeffs = n_bands * n_bands;
    constexpr int max_depth = 50;
    struct SHSample
    {
        Vector3f direction;
        float theta, phi;
        std::array<float, n_coeffs> Ylm;
    };
    typedef std::array<Vector3f, n_coeffs> SHCoefficients;
}

#endif