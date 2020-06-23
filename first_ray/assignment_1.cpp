#include "assignment_1.hpp"
#include "material.h"

Vector3f ass_1_renderer::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const double &prev_bsdf_pdf, sampler &random_sampler)
{
    hit_record hrec;
    auto &world = scene->world;
    auto &lights = scene->lights;
    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        Vector3f Le(0, 0, 0);
        scatter_record srec(hrec);

        if (hrec.mat_ptr->scatter(r, hrec, srec, random_sampler.get3d()))
        {
            /* Direct light sampling */
            for (auto *light : lights.list)
            {
                hit_record lrec;
                Vector3f offset_origin = hrec.p + (EPSILON * hrec.normal);
                Vector3f to_light = light->sample_direct(lrec, offset_origin, random_sampler.get2d());

                ray shadow_ray = ray(offset_origin, to_light);

                Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(shadow_ray, hrec, to_light);
                const double cos_wi = std::max(0.0, dot(hrec.normal, unit_vector(to_light)));

                const double light_pdf = light->pdf_direct_sampling(hrec, to_light);

                Le += surface_bsdf * lrec.mat_ptr->emitted(shadow_ray, lrec) * cos_wi / light_pdf;
            }
        }
        return Le;
    }

    return scene->env_map->eval(r, prev_hrec, depth);
}

void ass_1_renderer::Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf)
{
    tf.parallel_for(film_viewer->ny - 1, -1, -1, [=](int y)
        {
            static thread_local sampler random_sampler(y * 39);
            for (int x = 0; x < film_viewer->nx; x++)
            {
                if (film_viewer->to_exit) break;
                Vector3f col(0.0, 0.0, 0.0);
                for (int s = 0; s < film_viewer->ns; s++)
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
            std::cout << ".";
        });
}
