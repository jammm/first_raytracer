#ifndef PATH_PRT_H
#define PATH_PRT_H

#include "integrator.h"
#include "SHSample.h"
#include "material.h"

/* Path tracer using Precomputed Radiance Transfer for environment mapping with MIS using power heuristic */

struct path_prt
{
    // TODO
    // Convert from recursive to iterative
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
		const float &prev_bsdf_pdf);

	// Here, n_coeffs = n_bands*n_bands and n_samples = sqrt_n_samples*sqrt_n_samples
	void SHProject(int n_samples, int n_coeffs, const std::array<PRT::SHSample, PRT::n_coeffs> samples, double result[]);

	path_prt()
	{
		env_map = std::make_unique<environment_map>("");
	}

	std::unique_ptr<environment_map> env_map;

};

#endif