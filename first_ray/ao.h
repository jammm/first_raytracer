#ifndef AO_H
#define AO_H

#include "integrator.h"

/* Simple ambient occlusion integrator */

struct ao
{
    // TODO
    // Convert from recursive to iterative
    Vector3f Li(const ray& r, Scene* scene, const int& depth, const hit_record& prev_hrec,
        const float& prev_bsdf_pdf);
};


#endif