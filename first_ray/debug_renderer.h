#ifndef DEBUG_RENDERER_H
#define DEBUG_RENDERER_H

#include "integrator.h"

struct normals_renderer
{
    Vector3f Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
        const double &prev_bsdf_pdf, sampler &random_sampler)
	{
		hit_record hrec;
		if (scene->world->hit(r, EPSILON, FLT_MAX, hrec))
		{
			return hrec.normal;
		}
		return scene->env_map->eval(r, prev_hrec, depth);
	}

    void Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf)
    {
        tf.parallel_for(film_viewer->ny - 1, -1, -1, [=](int y)
            {
                static thread_local sampler random_sampler(y);
                for (uint_fast64_t x = 0; x < film_viewer->nx; x++)
                {
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
                }
                std::cout << ".";
            });
    }

	constexpr static bool using_custom_viewer = false;
};

#endif
