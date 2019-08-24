#include "path_prt.h"
#include "prt.h"
#include "triangle.h"
#include <taskflow/taskflow.hpp>

using namespace PRT;

path_prt::path_prt(Scene *scene, int &n_samples) : scene(scene), n_samples(n_samples)
{
    samples = PreComputeSamples(std::sqrt(n_samples), n_bands);
    SH_project_environment();
    SH_project_shadowed_diffuse_transfer();

    // Acual rendering after PRT only needs 1spp
    //n_samples = 1;
}

// Compute per-vertex coefficients for transfer function
// Transfer function encodes clamped cosine and diffuse BRDF terms together
void path_prt::SH_project_unshadowed_diffuse_transfer()
{
    hitable_list* world = dynamic_cast<hitable_list*>(scene->world.get());
    assert(world != nullptr);
    for (auto &i : world->list)
    {
        triangle* tri = dynamic_cast<triangle*>(i);
        if (tri == nullptr) continue;

        // Initialize per-vertex SH coefficients per triangle
        tri->coeffs.reset(new SHCoefficients[3]);
        for (int idx = 0; idx < 3; ++idx)
        {
            const Vector3f& v = tri->mesh->vertices[tri->V[idx]];
            const Vector3f& n = tri->mesh->normals[tri->V[idx]];
            const Point2f uv = tri->mesh->uv[tri->V[idx]];
            for (int i = 0; i < n_samples; ++i)
            {
                const Vector3f& direction = samples[i].direction;
                const float cosine = std::max(dot(n, direction), 0.0f);
                if (cosine == 0.0f) continue;
                hit_record rec;
                rec.p = v;
                rec.u = uv.x;
                rec.v = uv.y;
                const Vector3f albedo = tri->mat_ptr->get_albedo(rec);
                for (int n = 0; n < n_coeffs; ++n)
                {
                    const float value = cosine * samples[i].Ylm[n];
                    tri->coeffs[idx][n] += albedo * value / M_PI;
                }
            }
            // Divide the result by weight and number of samples
            const float factor = 4.0 * M_PI / n_samples;
            for (int i = 0; i < n_coeffs; ++i)
            {
                tri->coeffs[idx][i] *= factor;
            }
        }
    }
}

// Compute per-vertex coefficients for shadowed transfer function
void path_prt::SH_project_shadowed_diffuse_transfer()
{
    hitable_list* world = dynamic_cast<hitable_list*>(scene->world.get());
    assert(world != nullptr);
    tf::Taskflow tf;
    int world_size = world->list.size();
    tf.parallel_for(0, world_size, 1, [&](int obj_idx)
    {
        triangle* tri = dynamic_cast<triangle*>(world->list[obj_idx]);
        if (tri != nullptr)
        {
            // Initialize per-vertex SH coefficients per triangle
            tri->coeffs.reset(new SHCoefficients[3]);
            for (int idx = 0; idx < 3; ++idx)
            {
                const Vector3f& v = tri->mesh->vertices[tri->V[idx]];
                const Vector3f& n = tri->mesh->normals[tri->V[idx]];
                const Point2f uv = tri->mesh->uv[tri->V[idx]];
                for (int i = 0; i < n_samples; ++i)
                {
                    const Vector3f& direction = samples[i].direction;
                    const float cosine = std::max(dot(n, direction), 0.0f);
                    if (cosine == 0.0f) continue;
                    const ray r(v + (EPSILON * n), direction);
                    hit_record rec;
                    if (!scene->world->hit(r, EPSILON, FLT_MAX, rec))
                    {
                        rec.p = v;
                        rec.u = uv.x;
                        rec.v = uv.y;
                        const Vector3f albedo = tri->mat_ptr->get_albedo(rec);
                        for (int n = 0; n < n_coeffs; ++n)
                        {
                            const float value = cosine * samples[i].Ylm[n];
                            tri->coeffs[idx][n] += albedo * value / M_PI;
                        }
                    }
                }
                // Divide the result by weight and number of samples
                const float factor = 4.0 * M_PI / n_samples;
                for (int i = 0; i < n_coeffs; ++i)
                {
                    tri->coeffs[idx][i] *= factor;
                    if (!coeffs_buffer.empty())
                    {
                        coeffs_buffer[obj_idx + idx][0][i] = tri->coeffs[idx][i];
                    }
                }
            }
        }
    });
    tf.dispatch().get();
}

// SH projection for diffuse inter-reflection
// It calculates coefficients for each vertex, adding
// coefficients calculated from direct lighting (diffuse shadowed SH projection).
// The process is repeated again for the next bounce by storing calcualted coefficients
// and using them as the base for the next bounce.
// Finally, every vertex's coefficients are summed together from all bounces,
// giving the final transfer vector.

void path_prt::SH_project_full_global_illumination()
{
    hitable_list* world = dynamic_cast<hitable_list*>(scene->world.get());
    assert(world != nullptr);
    // Allocate memory for coefficients buffer
    // coeffs_buffer contains coefficients for each vertex stored separately per bounce
    coeffs_buffer = std::vector<SHCoefficients[max_depth]>(world->list.size() * 3);
    // Fill coefficients on depth 0 (direct illumination) using shadowed SH projection
    SH_project_shadowed_diffuse_transfer();
    const unsigned int world_size = world->list.size();

    tf::Taskflow tf;
    // For upto 10 bounces
    // TODO: Use max_depth instead of hardcoding 10
    tf.parallel_for(0, 10, 1, [&](unsigned int depth)
    {
        std::cout << depth << std::endl;
        // Loop through all triangles in scene
        for (unsigned int obj_idx = 0; obj_idx < world_size; ++obj_idx)
        {
            // Get current object. Ignore if object isn't triangle
            triangle* tri = dynamic_cast<triangle*>(world->list[obj_idx]);

            if (tri == nullptr) continue;
            // Process each vertex of triangle
            for (int idx = 0; idx < 3; ++idx)
            {
                const Vector3f& v = tri->mesh->vertices[tri->V[idx]];
                const Vector3f& n = tri->mesh->normals[tri->V[idx]];
                const Point2f uv = tri->mesh->uv[tri->V[idx]];
                for (int i = 0; i < n_samples; ++i)
                {
                    // Get precomputed sample direction
                    const Vector3f& direction = samples[i].direction;
                    const ray r(v + (EPSILON * n), direction);

                    hit_record rec;
                    // Cast a ray to the scene
                    if (scene->world->hit(r, EPSILON, FLT_MAX, rec))
                    {
                        const float cosine = std::max(dot(n, direction), 0.0f);
                        if (cosine == 0.0f) continue;
                        // Ray is within hemisphere
                        auto& tri_coeffs = dynamic_cast<triangle*>(rec.obj)->coeffs;
                        SHCoefficients cur_coeffs = {};

                        // Interpolate coefficients on hit point
                        for (int n = 0; n < n_coeffs; ++n)
                        {
                            // Use coefficients from previous bounce
                            cur_coeffs[n] += (1 - uv.x - uv.y) * coeffs_buffer[obj_idx][depth-1][n]
                                                        + uv.x * coeffs_buffer[obj_idx + 1][depth-1][n]
                                                        + uv.y * coeffs_buffer[obj_idx + 2][depth-1][n];
                        }
                        const Point2f &uv = rec.uv;
                        rec.p = v;
                        rec.u = uv.x;
                        rec.v = uv.y;
                        const Vector3f albedo = tri->mat_ptr->get_albedo(rec);
                        // Calculate coefficients for current bounce
                        for (int n = 0; n < n_coeffs; ++n)
                        {
                            const Vector3f value = cosine * cur_coeffs[n];
                            coeffs_buffer[obj_idx + idx][depth][n] += albedo * value / M_PI;
                        }
                    }
                }

            }
        }
        // Normalize coefficients
        for (unsigned int obj_idx = 0; obj_idx < world_size; ++obj_idx)
        {
            // Divide the result by weight and number of samples
            const float factor = 4.0 * M_PI / n_samples;
            for (int n = 0; n < n_coeffs; ++n)
            {
                coeffs_buffer[obj_idx][depth][n] *= factor;
                coeffs_buffer[obj_idx + 1][depth][n] *= factor;
                coeffs_buffer[obj_idx + 2][depth][n] *= factor;
            }
        }
    });
    tf.dispatch().get();

    // Sum all coefficients for current bounce
    for (unsigned int obj_idx = 0; obj_idx < world_size; ++obj_idx)
    {
        triangle* tri = dynamic_cast<triangle*>(world->list[obj_idx]);
        if (tri == nullptr) continue;

        // Each triangle vertex has its own coefficient vector
        for (unsigned int depth = 1; depth <= max_depth; ++depth)
            for (int idx = 0; idx < 3; ++idx)
                for (int n = 0; n < n_coeffs; ++n)
                    tri->coeffs[idx][n] += coeffs_buffer[obj_idx + idx][depth][n];
    }
}

// Here, n_coeffs = n_bands*n_bands and n_samples = sqrt_n_samples*sqrt_n_samples
void path_prt::SH_project_environment()
{
    Li_coeffs = {};
    // For each sample
    for (int i = 0; i < n_samples; ++i)
    {
        hit_record envmap_rec;
        ray r(Vector3f(0, 0, 0), samples[i].direction);
        for (int n = 0; n < n_coeffs; ++n)
        {
            Li_coeffs[n] += scene->env_map->eval(r, envmap_rec, -1) * samples[i].Ylm[n];
        }
    }
    // Divide the result by weight and number of samples
    const float factor = 4.0 * M_PI / n_samples;
    for (int i = 0; i < n_coeffs; ++i)
    {
        Li_coeffs[i] *= factor;
    }
}

// TODO
// Convert from recursive to iterative
Vector3f path_prt::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const float &prev_bsdf_pdf)
{
    hit_record hrec;
    auto &world = scene->world;
    auto &lights = scene->lights;
    Vector3f result(0.0f, 0.0f, 0.0f);

    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec;
        if (hrec.mat_ptr->scatter(r, hrec, srec))
        {
            auto &tri_coeffs = dynamic_cast<triangle*>(hrec.obj)->coeffs;
            Point2f &uv = hrec.uv;
            SHCoefficients tr_coeffs = {};

            for (int i = 0; i < n_coeffs; ++i)
            {
                tr_coeffs[i] += (1 - uv.x - uv.y) * tri_coeffs[0][i] + uv.x * tri_coeffs[1][i] + uv.y * tri_coeffs[2][i];
                result += tr_coeffs[i] * Li_coeffs[i];
            }
            return result;
        }
    }
    return scene->env_map->eval(r, hrec, depth);
}

