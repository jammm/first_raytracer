#ifndef DEBUG_RENDERER_H
#define DEBUG_RENDERER_H

#include "integrator.h"

struct normals_renderer
{
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
        const float &prev_bsdf_pdf, sampler &random_sampler)
	{
		hit_record hrec;
		if (scene->world->hit(r, EPSILON, FLT_MAX, hrec))
		{
			return hrec.normal;
		}
		return scene->env_map->eval(r, prev_hrec, depth);
	}
};

#endif