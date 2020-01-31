#ifndef PATH_H
#define PATH_H

#include "integrator.h"

/* Path tracer with MIS using power heuristic */

struct path
{
    // TODO
    // Convert from recursive to iterative
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
        const float &prev_bsdf_pdf, sampler &random_sampler);

    void Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf);
};

#endif