#ifndef PRT_H
#define PRT_H

#include "util.h"
#include "geometry.h"
#include "SHSample.h"


namespace PRT
{
    constexpr unsigned int factorial(const unsigned int &n)
    {
        if (n < 2) return 1;
        return n * factorial(n - 1);
    }

	// Renormalisation constant for SH function
	double K(int l, int m) {
		double temp = ((2.0 * l + 1.0) * factorial(l - m)) / (4.0 * M_PI * factorial(l + m));   // Here, you can use a precomputed table for factorials
		return sqrt(temp);
	}

	// Evaluate an Associated Legendre Polynomial P(l,m,x) at x
	// For more, see “Numerical Methods in C: The Art of Scientific Computing”, Cambridge University Press, 1992, pp 252-254 
	double P(int l, int m, double x) {
		double pmm = 1.0;
		if (m > 0) {
			double somx2 = sqrt((1.0 - x) * (1.0 + x));
			double fact = 1.0;
			for (int i = 1; i <= m; i++) {
				pmm *= (-fact) * somx2;
				fact += 2.0;
			}
		}
		if (l == m)
			return pmm;

		double pmmp1 = x * (2.0 * m + 1.0) * pmm;
		if (l == m + 1)
			return pmmp1;

		double pll = 0.0;
		for (int ll = m + 2; ll <= l; ++ll) {
			pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
			pmm = pmmp1;
			pmmp1 = pll;
		}

		return pll;
	}

	// Returns a point sample of a Spherical Harmonic basis function
	// l is the band, range [0..N]
	// m in the range [-l..l]
	// theta in the range [0..Pi]
	// phi in the range [0..2*Pi]
	double EstimateSH(int l, int m, double theta, double phi) {
		const double sqrt2 = sqrt(2.0);
		if (m == 0)        return K(l, 0) * P(l, m, cos(theta));
		else if (m > 0)    return sqrt2 * K(l, m) * cos(m * phi) * P(l, m, cos(theta));
		else                return sqrt2 * K(l, -m) * sin(-m * phi) * P(l, -m, cos(theta));
	}

	// Fills an N*N*2 array with uniformly distributed samples across the sphere using jittered stratification
	std::unique_ptr<SHSample[]> PreComputeSamples(int sqrt_n_samples, int n_bands) {
		double oneoverN = 1.0 / sqrt_n_samples;

		// Create samples array
		auto samples = std::make_unique<SHSample[]>(sqrt_n_samples * sqrt_n_samples);
		sampler random_sampler;

		for (int a = 0, i=0; a < sqrt_n_samples; a++) {
			for (int b = 0; b < sqrt_n_samples; b++, i++) {
				// Generate unbiased distribution of spherical coords
				double x = (a + random_sampler.get1d()) * oneoverN;           // Do not reuse results
				double y = (b + random_sampler.get1d()) * oneoverN;           // Each sample must be random!
				double theta = 2.0 * acos(sqrt(1.0 - x));   // Uniform sampling on theta
				double phi = 2.0 * M_PI * y;                      // Uniform sampling on phi

				// Convert spherical coords to unit vector
				samples[i].direction = Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
				samples[i].theta = theta;
				samples[i].phi = phi;

				// precompute all SH coefficients for this sample
				for (int l = 0; l < n_bands; ++l) {
					for (int m = -l; m <= l; ++m) {
						int index = l * (l + 1) + m;
						samples[i].Ylm[index] = EstimateSH(l, m, theta, phi);
					}
				}
			}
		}
		return samples;
	}
}



#endif
