#include "ao.h"
#include "material.h"

Vector3f ao::Li(const ray& r, Scene* scene, const int& depth, const hit_record& prev_hrec,
    const float& prev_bsdf_pdf)
{
    hit_record hrec;
    auto& world = scene->world;
    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec;

        if (hrec.mat_ptr->scatter(r, hrec, srec))
        {
            /* Sample BSDF to generate next ray direction for visibility check */
            hrec.p = hrec.p;
            ray shadow_ray(hrec.p, srec.pdf_ptr->generate());
            // Find approximate max. radius from bounding box
            // TODO: evaluate bounding sphere for every hitable geometry
            aabb bbox;
            scene->world->bounding_box(0.0f, 0.0f, bbox);
            const float t_max = bbox.size.y() * 0.50f;
            if (world->hit(shadow_ray, EPSILON, t_max, hrec))
                return Vector3f(0.0f, 0.0f, 0.0f);
        }
    }
    return scene->env_map->eval(r, prev_hrec, depth);
}