#ifndef PATH_PRT_H
#define PATH_PRT_H

#include "integrator.h"
#include "SHSample.h"
#include "material.h"

/* Path tracer using Precomputed Radiance Transfer for environment mapping with MIS using power heuristic */

struct path_prt
{
	path_prt(Scene *scene, int &n_samples);
    
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
		const float &prev_bsdf_pdf);

	void SH_project_environment();
	void SH_project_unshadowed_diffuse_transfer();
	void SH_project_shadowed_diffuse_transfer();

	const Scene *scene;
	int n_samples;
	std::unique_ptr<Vector3f[]> Li_coeffs;
	std::unique_ptr<PRT::SHSample[]> samples;
};

#endif