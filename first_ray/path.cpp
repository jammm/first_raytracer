#include "path.h"
#include "material.h"

Vector3f path::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const double &prev_bsdf_pdf, sampler &random_sampler)
{
    hit_record hrec;
    auto &world = scene->world;
    auto &lights = scene->lights;
    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec(hrec);
        Vector3f Le = hrec.mat_ptr->emitted(r, hrec);

        /* If we hit a light source, weight its contribution */
        if (((Le.r() != 0.0) || (Le.g() != 0.0) || (Le.b() != 0.0)))
        {
            if ((depth == 0)
                || (dynamic_cast<modified_phong*>(prev_hrec.mat_ptr) != nullptr)
                || (dynamic_cast<metal*>(prev_hrec.mat_ptr) != nullptr)
                || (dynamic_cast<dielectric*>(prev_hrec.mat_ptr) != nullptr))
                return Le;
            // Start with checking if camera ray hits a light source
            const double cos_wo = dot(hrec.normal, -unit_vector(r.direction()));
            double distance_squared = (hrec.p - prev_hrec.p).squared_length();
            if (distance_squared <= EPSILON) distance_squared = EPSILON;

            double surface_bsdf_pdf = prev_bsdf_pdf;
            const double light_pdf = hrec.obj->pdf_direct_sampling(hrec, r.direction()) * distance_squared / abs(cos_wo);

            const double weight = miWeight(surface_bsdf_pdf, light_pdf);

            return Le * weight;
        }

        if (depth <= 33 && hrec.mat_ptr->scatter(r, hrec, srec, random_sampler.get3d()))
        {
            /* Direct light sampling */
            const int index = lights.pick_sample(random_sampler.get1d());
            if (index >= 0 && (dynamic_cast<dielectric*>(hrec.mat_ptr) == nullptr))
            {
                /* Sample a random light source */
                hit_record lrec;
                Vector3f offset_origin = hrec.p + (EPSILON * hrec.normal);
                Vector3f to_light = lights[index]->sample_direct(lrec, offset_origin, random_sampler.get2d());
                const double dist_to_light = to_light.length();

                ray shadow_ray = ray(offset_origin, to_light);

                if (!world->hit(shadow_ray, EPSILON, 1 - SHADOW_EPSILON, lrec))
                {
                    to_light.make_unit_vector();
                    shadow_ray.d = to_light;
                    Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(shadow_ray, hrec, to_light);
                    // Calculate geometry term
                    const double cos_wi = dot(hrec.normal, unit_vector(to_light));
                    const double cos_wo = dot(lrec.normal, -unit_vector(to_light));
                    if (cos_wo != 0)
                    {
                        double distance_squared = dist_to_light * dist_to_light;

                        if (!srec.is_specular)
                            surface_bsdf *= cos_wi;

                        const double light_pdf = lights[index]->pdf_direct_sampling(hrec, to_light) * distance_squared / abs(cos_wo);
                        // Visibility term is always 1
                        // because of the invariant imposed on these objects by the if above.
                        const double surface_bsdf_pdf = srec.pdf_ptr->value(hrec, to_light);

                        if (distance_squared <= EPSILON) distance_squared = EPSILON;

                        const double weight = miWeight(light_pdf, surface_bsdf_pdf);

                        Le += lrec.mat_ptr->emitted(shadow_ray, lrec) * surface_bsdf * weight / light_pdf;
                    }
                }
            }
            if (srec.is_specular)
            {
                double surface_bsdf_pdf = srec.pdf_ptr ? srec.pdf_ptr->value(hrec, srec.specular_ray.direction()) : 1.0;
                if (srec.sampled_pdf > 0.0) surface_bsdf_pdf = srec.sampled_pdf;
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(r, hrec, srec.specular_ray.direction());
                assert(std::isfinite(surface_bsdf[0])
                    && std::isfinite(surface_bsdf[1])
                    && std::isfinite(surface_bsdf[2]));
                assert(std::isfinite(surface_bsdf_pdf));
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                //const double cos_wo = abs(dot(hrec.normal, unit_vector(srec.specular_ray.direction())));
                const bool outside = dot(hrec.normal, srec.specular_ray.d) > 0;
                srec.specular_ray.o = outside ? (srec.specular_ray.o + (EPSILON * hrec.normal)) : (srec.specular_ray.o - (EPSILON * hrec.normal));
                return Le + surface_bsdf * Li(srec.specular_ray, scene, depth + 1, hrec, surface_bsdf_pdf, random_sampler) / surface_bsdf_pdf;
            }
            else
            {
                /* Sample BSDF to generate next ray direction for indirect lighting */
                ray wo(hrec.p + EPSILON * (hrec.normal), srec.pdf_ptr->generate(random_sampler.get2d(), srec));
                const double surface_bsdf_pdf = srec.pdf_ptr->value(hrec, wo.direction());
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(wo, hrec, wo.direction());
                /* Reject current path in case the ray is on the wrong side of the surface (BRDF is 0 as ray is pointing away from the hemisphere )*/
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                const double cos_wo = abs(dot(hrec.normal, unit_vector(wo.direction())));

                return Le + surface_bsdf * Li(wo, scene, depth + 1, hrec, surface_bsdf_pdf, random_sampler) * cos_wo / surface_bsdf_pdf;
            }
        }
        return Le;
    }

    return scene->env_map->eval(r, prev_hrec, depth);
}

void path::Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf)
{
    tf.parallel_for(210 - 1, 100, -1, [=](int y)
        {
            static thread_local sampler random_sampler(y * 39);
            for (int x = 0; x < film_viewer->nx; x++)
            {
                //if (x == 360 && y == 127)
                {
                    if (film_viewer->to_exit) break;
                    Vector3f col(0.0, 0.0, 0.0);
                    for (uint_fast64_t s = 0; s < film_viewer->ns; s++)
                    {
                        double u = double(x + random_sampler.get1d()) / double(film_viewer->nx);
                        double v = double(y + random_sampler.get1d()) / double(film_viewer->ny);
                        ray r = scene->cam.get_ray(u, v, random_sampler.get2d());
                        hit_record hrec;
                        // Compute a sample
                        const Vector3f sample = Li(r, scene, 0, hrec, 0.0, random_sampler);
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
