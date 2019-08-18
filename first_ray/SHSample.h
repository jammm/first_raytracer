#ifndef SH_SAMPLE_H
#define SH_SAMPLE_H

#include "geometry.h"
#include "util.h"

namespace PRT
{

	constexpr int n_coeffs = 3;
	struct SHSample
	{
		Vector3f  direction;
		float   theta, phi;
		std::array<float, n_coeffs> Ylm;
	};
}

#endif