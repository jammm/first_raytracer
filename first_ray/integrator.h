#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "hitable_list.h"
#include "geometry.h"
#include "Scene.h"

template <typename integrator>
struct renderer : public integrator
{
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
        const float &prev_bsdf_pdf)
    {
        return integrator::Li(r, scene, depth, prev_hrec, prev_bsdf_pdf);
    };
};

#endif
