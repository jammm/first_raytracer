#ifndef PRT_H
#define PRT_H

#include "util.h"
#include "geometry.h"

namespace PRT
{
    constexpr int n_bands = 3;
    struct SHSample {
        Vector3f sph;
        Vector3f vec;
        double *coeff;
    };

    double P(int l, int m, double x)
    {
        // evaluate an Associated Legendre Polynomial P(l,m,x) at x
        double pmm = 1.0;
        if (m > 0) {
            double somx2 = sqrt((1.0 - x)*(1.0 + x));
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

    double K(int l, int m)
    {
        // renormalisation constant for SH function
        double temp = ((2.0*l + 1.0)*factorial(l - m)) / (4.0*M_PI*factorial(l + m));
        return sqrt(temp);
    }
    double SH(int l, int m, double theta, double phi)
    {
        // return a point sample of a Spherical Harmonic basis function
        // l is the band, range [0..N]
        // m in the range [-l..l]
        // theta in the range [0..Pi]
        // phi in the range [0..2*Pi]
        const double sqrt2 = sqrt(2.0);
        if (m == 0) return K(l, 0)*P(l, m, cos(theta));
        else if (m > 0) return sqrt2 * K(l, m)*cos(m*phi)*P(l, m, cos(theta));
        else return sqrt2 * K(l, -m)*sin(-m * phi)*P(l, -m, cos(theta));
    }
    void SH_setup_spherical_samples(std::vector<SHSample> samples, int sqrt_n_samples)
    {
        // fill an N*N*2 array with uniformly distributed
        // samples across the sphere using jittered stratification
        int i = 0; // array index
        double oneoverN = 1.0 / sqrt_n_samples;
        for (int a = 0; a < sqrt_n_samples; a++) {
            for (int b = 0; b < sqrt_n_samples; b++) {
                // generate unbiased distribution of spherical coords
                const double x = (a + gen_cano_rand()) * oneoverN; // do not reuse results
                const double y = (b + gen_cano_rand()) * oneoverN; // each sample must be random
                const double theta = 2.0 * acos(sqrt(1.0 - x));
                const double phi = 2.0 * M_PI * y;
                samples[i].sph = Vector3f(theta, phi, 1.0);
                const double sin_theta = sin(theta);
                // convert spherical coords to unit vector
                Vector3f vec(sin_theta*cos(phi), sin_theta*sin(phi), cos(theta));
                samples[i].vec = vec;
                // precompute all SH coefficients for this sample
                for (int l = 0; l < n_bands; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        int index = l * (l + 1) + m;
                        samples[i].coeff[index] = SH(l, m, theta, phi);
                    }
                }
                ++i;
            }
        }
    }
}



#endif;