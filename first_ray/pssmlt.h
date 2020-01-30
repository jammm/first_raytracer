#ifndef PATH_H
#define PATH_H

#include "integrator.h"

/* Implements the PSSMLT paper "Simple and Robust Mutation Strategy for Metropolis Light Transport Algorithm"
   - Csaba Kelemen and László Szirmay-Kalos */

struct pssmlt
{
    // TODO
    // Convert from recursive to iterative
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
<<<<<<< HEAD
        const float &prev_bsdf_pdf, sampler &random_sampler);

    void Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf);
=======
        const float &prev_bsdf_pdf);
>>>>>>> 789dbd6... Refactor viewer into a separate class
};

#endif