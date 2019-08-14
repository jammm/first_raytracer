#ifndef PRT_H
#define PRT_H

#include "util.h"
#include "geometry.h"

#include <gcem/gcem.hpp>
#include <array>

namespace PRT
{

    template<unsigned int n_bands = 3>
    struct SHSample {
        Vector3f sph;
        Vector3f vec;

        constexpr SHSample() = default;

        std::array<double, n_bands> coeff;
    };

    constexpr unsigned int factorial(const unsigned int &n)
    {
        if (n < 2) return 1;
        return n * factorial(n - 1);
    }

    constexpr double P(int l, int m, double x)
    {
        // evaluate an Associated Legendre Polynomial P(l,m,x) at x
        double pmm = 1.0;
        if (m > 0) {
            double somx2 = gcem::sqrt((1.0 - x)*(1.0 + x));
            double fact = 1.0;
            for (int i = 1; i <= m; i++) {
                pmm *= (-fact) * somx2;
                fact += 2.0;
            }
        }
        if (l == m) return pmm;
        double pmmp1 = x * (2.0*m + 1.0) * pmm;
        if (l == m + 1) return pmmp1;
        double pll = 0.0;
        for (int ll = m + 2; ll <= l; ++ll) {
            pll = ((2.0*ll - 1.0)*x*pmmp1 - (ll + m - 1.0)*pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }
        return pll;
    }

    constexpr double K(const int &l, const int m)
    {
        // renormalisation constant for SH function
        const double temp = ((2.0*l + 1.0)*factorial(l - m)) / (4.0*M_PI*factorial(l + m));
        return gcem::sqrt(temp);
    }

    constexpr double SH(int l, int m, double theta, double phi)
    {
        // return a point sample of a Spherical Harmonic basis function
        // l is the band, range [0..N]
        // m in the range [-l..l]
        // theta in the range [0..Pi]
        // phi in the range [0..2*Pi]
        const double sqrt2 = gcem::sqrt(2.0);
        if (m == 0) return K(l, 0)*P(l, m, gcem::cos(theta));
        else if (m > 0) return sqrt2 * K(l, m)*gcem::cos(m*phi)*P(l, m, gcem::cos(theta));
        else return sqrt2 * K(l, -m)*gcem::sin(-m * phi)*P(l, -m, gcem::cos(theta));
    }

    template<unsigned int n_bands, unsigned int n_samples>
    constexpr auto SH_setup_spherical_samples() -> decltype(auto)
    {
        std::array<SHSample<n_bands>, n_samples> samples;
        const auto sqrt_n_samples = gcem::sqrt(n_samples);
        // fill an N*N*2 array with uniformly distributed
        // samples across the sphere using jittered stratification
        double oneoverN = 1.0 / sqrt_n_samples;
        for (int a = 0, i = 0; a < sqrt_n_samples; a++) {
            for (int b = 0; b < sqrt_n_samples; b++, i++) {
                // generate unbiased distribution of spherical coords
                const double x = (a + 0) * oneoverN; // do not reuse results
                const double y = (b + 0) * oneoverN; // each sample must be random
                const double theta = 2.0 * gcem::acos(sqrt(1.0 - x));
                const double phi = 2.0 * M_PI * y;
                samples[i].sph = Vector3f(theta, phi, 1.0);
                const double sin_theta = gcem::sin(theta);
                // convert spherical coords to unit vector
                Vector3f vec(sin_theta*gcem::cos(phi), sin_theta*gcem::sin(phi), gcem::cos(theta));
                samples[i].vec = vec;
                // precompute all SH coefficients for this sample
                for (int l = 0; l < n_bands; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        int index = l * (l + 1) + m;
                        samples[i].coeff[index] = SH(l, m, theta, phi);
                    }
                }
            }
        }

        return samples;
    }

    constexpr auto __data = SH_setup_spherical_samples<3, 10000>();
}



#endif;