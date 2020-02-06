#include "path.h"
#include "material.h"

Vector3f path::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const float &prev_bsdf_pdf, sampler &random_sampler)
{
    hit_record hrec;
    auto &world = scene->world;
    auto &lights = scene->lights;
    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec;
        Vector3f Le = hrec.mat_ptr->emitted(r, hrec);

        /* If we hit a light source, weight its contribution */
        if (((Le.r() != 0.0f) || (Le.g() != 0.0f) || (Le.b() != 0.0f)))
        {
            if ((depth == 0)
                || (dynamic_cast<modified_phong*>(prev_hrec.mat_ptr) != nullptr)
                || (dynamic_cast<metal*>(prev_hrec.mat_ptr) != nullptr))
                return Le;
            // Start with checking if camera ray hits a light source
            const float cos_wo = std::max(dot(hrec.normal, -unit_vector(r.direction())), 0.0f);
            float distance_squared = (hrec.p - prev_hrec.p).squared_length();
            if (distance_squared <= EPSILON) distance_squared = EPSILON;

            const float surface_bsdf_pdf = prev_bsdf_pdf * cos_wo / distance_squared;
            const float light_pdf = lights.pdf_direct_sampling(hrec, r.direction());

            const float weight = miWeight(surface_bsdf_pdf, light_pdf);

            return Le * weight;
        }

        if (depth <= 50 && hrec.mat_ptr->scatter(r, hrec, srec, random_sampler.get3d()))
        {
            if (srec.is_specular)
            {
                const float surface_bsdf_pdf = srec.pdf_ptr ? srec.pdf_ptr->value(hrec, srec.specular_ray.direction()) : 1.0f;
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(r, hrec, srec.specular_ray.direction());
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                //const float cos_wi = abs(dot(hrec.normal, unit_vector(srec.specular_ray.direction())));
                return surface_bsdf * Li(srec.specular_ray, scene, depth + 1, hrec, surface_bsdf_pdf, random_sampler) / surface_bsdf_pdf;
            }
            else
            {
                /* Direct light sampling */
                const int index = lights.pick_sample(random_sampler.get1d());
                if (index >= 0)
                {
                    /* Sample a random light source */
                    hit_record lrec;
                    Vector3f offset_origin = hrec.p + (EPSILON * hrec.normal);
                    Vector3f to_light = lights[index]->sample_direct(lrec, offset_origin, random_sampler.get2d());
                    const float dist_to_light = to_light.length();
                    //to_light.make_unit_vector();

                    ray shadow_ray = ray(offset_origin, to_light);

                    if (!world->hit(shadow_ray, EPSILON, 1 - SHADOW_EPSILON, lrec))
                    {
                        const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(shadow_ray, hrec, to_light);
                        // Calculate geometry term
                        const float cos_wi = std::abs(dot(hrec.normal, unit_vector(to_light)));
                        const float cos_wo = std::max(dot(lrec.normal, -unit_vector(to_light)), 0.0f);
                        if (cos_wo != 0)
                        {
                            float distance_squared = dist_to_light * dist_to_light;

                            const float light_pdf = lights.pdf_direct_sampling(hrec, to_light);
                            // Visibility term is always 1
                            // because of the invariant imposed on these objects by the if above.
                            const float surface_bsdf_pdf = srec.pdf_ptr->value(hrec, to_light) * cos_wo / distance_squared;

                            if (distance_squared <= EPSILON) distance_squared = EPSILON;

                            const float G = cos_wi * cos_wo / distance_squared;

                            const float weight = miWeight(light_pdf, surface_bsdf_pdf);

                            Le += lights.list_size * lrec.mat_ptr->emitted(shadow_ray, lrec) * surface_bsdf * G * weight / light_pdf;
                        }
                    }
                }
                /* Sample BSDF to generate next ray direction for indirect lighting */
                hrec.p = hrec.p + (EPSILON * hrec.normal);
                ray wo(hrec.p, srec.pdf_ptr->generate(random_sampler.get2d(), hrec));
                const float surface_bsdf_pdf = srec.pdf_ptr->value(hrec, wo.direction());
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(wo, hrec, wo.direction());
                /* Reject current path in case the ray is on the wrong side of the surface (BRDF is 0 as ray is pointing away from the hemisphere )*/
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                const float cos_wo = abs(dot(hrec.normal, unit_vector(wo.direction())));

                return Le + surface_bsdf * Li(wo, scene, depth + 1, hrec, surface_bsdf_pdf, random_sampler) * cos_wo / surface_bsdf_pdf;
            }
        }
        return Le;
    }

    //Vector3f unit_direction = unit_vector(r.direction());
    //float t = 0.5*(unit_direction.y() + 1.0);
    //return (1.0 - t)*Vector3f(1.0, 1.0, 1.0) + t * Vector3f(0.5, 0.7, 1.0);
    return scene->env_map->eval(r, prev_hrec, depth);
}

void path::Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf)
{
    tf.parallel_for(film_viewer->ny - 1, 0, -1, [=](int y)
        {
            static thread_local sampler random_sampler(y);
            for (int x = 0; x < film_viewer->nx; x++)
            {
                //if (i <= 512)
                {
                    if (film_viewer->to_exit) break;
                    Vector3f col(0.0f, 0.0f, 0.0f);
                    for (int s = 0; s < film_viewer->ns; s++)
                    {
                        float u = float(x + random_sampler.get1d()) / float(film_viewer->nx);
                        float v = float(y + random_sampler.get1d()) / float(film_viewer->ny);
                        ray r = scene->cam.get_ray(u, v, random_sampler.get2d());
                        hit_record hrec;
                        // Compute a sample
                        const Vector3f sample = Li(r, scene, 0, hrec, 0.0f, random_sampler);
                        assert(std::isfinite(sample[0])
                            && std::isfinite(sample[1])
                            && std::isfinite(sample[2]));
                        col += sample;
                    }
                    // Splat this sample to film
                    film_viewer->add_sample(Vector2i(x, y), col);
                }
            }
            std::cout << ".";
        });
}