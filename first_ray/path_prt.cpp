#include "path_prt.h"
#include "prt.h"
#ifdef USE_SSE
#include "triangle_sse.hpp"
#else
#include "triangle.h"
#endif
#include <taskflow/taskflow.hpp>

using namespace PRT;

path_prt::path_prt(Scene *scene, int &n_samples) : scene(scene), n_samples(n_samples)
{
    samples = PreComputeSamples(std::sqrt(n_samples), n_bands);
    SH_project_environment();
    SH_project_full_global_illumination();

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
            const Vector2f uv = tri->mesh->uv[tri->V[idx]];
            for (int i = 0; i < n_samples; ++i)
            {
                const Vector3f& direction = samples[i].direction;
                const double cosine = std::max(dot(n, direction), 0.0);
                if (cosine == 0.0) continue;
                hit_record rec;
                rec.p = v;
                rec.u = uv.x;
                rec.v = uv.y;
                const Vector3f albedo = tri->mat_ptr->get_albedo(rec);
                for (int n = 0; n < n_coeffs; ++n)
                {
                    const double value = cosine * samples[i].Ylm[n];
                    tri->coeffs[idx][n] += albedo * value / M_PI;
                }
            }
            // Divide the result by weight and number of samples
            const double factor = 4.0 * M_PI / n_samples;
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
                const Vector2f uv = tri->mesh->uv[tri->V[idx]];
                for (int i = 0; i < n_samples; ++i)
                {
                    const Vector3f& direction = samples[i].direction;
                    const double cosine = std::max(dot(n, direction), 0.0);
                    if (cosine == 0.0) continue;
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
                            const double value = cosine * samples[i].Ylm[n];
                            tri->coeffs[idx][n] += albedo * value / M_PI;
                        }
                    }
                }
                // Divide the result by weight and number of samples
                const double factor = 4.0 * M_PI / n_samples;
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
    const unsigned int world_size = world->list.size();

    tf::Taskflow tf;
    // Loop through all triangles in scene
    tf.parallel_for(0U, world_size, 1U, [&](int obj_idx)
    {
        sampler random_sampler(obj_idx);
        // Get current object. Ignore if object isn't triangle
        triangle* tri = dynamic_cast<triangle*>(world->list[obj_idx]);

        //if (tri == nullptr) continue;
        // Process each vertex of triangle
        // Initialize per-vertex SH coefficients per triangle
        tri->coeffs.reset(new SHCoefficients[3]);
        for (int idx = 0; idx < 3; ++idx)
        {
            const Vector3f& v = tri->mesh->vertices[tri->V[idx]];
            const Vector3f n = tri->mesh->normals[tri->V[idx]];
            Vector3f cur_n = n;
            const Vector2f uv = tri->mesh->uv[tri->V[idx]];
            //Vector3f result(0.0, 0.0, 0.0);
            for (int i = 0; i < n_samples; ++i)
            {
                // Get precomputed sample direction
                Vector3f throughput(1.0, 1.0, 1.0);
                Vector3f result(0.0, 0.0, 0.0);
                hit_record hrec;
                hrec.p = v;
                hrec.normal = n;
                hrec.u = uv.x;
                hrec.v = uv.y;
                scatter_record srec(hrec);
                const Vector3f v_direction = cosine_pdf(n).generate(random_sampler.get2d(), srec);
                ray r(v + (EPSILON * n), v_direction);
                const Vector3f v_bsdf = tri->mat_ptr->eval_bsdf(r, hrec, r.direction());
                // Set initial BSDF to be used in case current vertex is directly visible from envmap
                Vector3f bsdf = v_bsdf;

                for (unsigned int depth = 0; depth < max_depth; ++depth)
                {
                    scatter_record srec(hrec);
                    // Cast a ray to the scene
                    if (scene->world->hit(r, EPSILON, FLT_MAX, hrec))
                    {
                        if (hrec.mat_ptr->scatter(r, hrec, srec, random_sampler.get3d()))
                        {
                            const double cos_wi = std::max(dot(hrec.normal, unit_vector(r.direction())), 0.0);
                            if (cos_wi == 0.0) continue;

                            /* Sample BSDF to generate next ray direction for indirect lighting */
                            hrec.p = hrec.p + (EPSILON * hrec.normal);
                            r = ray(hrec.p, srec.pdf_ptr->generate(random_sampler.get2d(), srec));
                            const double surface_bsdf_pdf = srec.pdf_ptr->value(hrec, r.direction());
                            bsdf = hrec.mat_ptr->eval_bsdf(r, hrec, r.direction());
                            /* Reject current path in case the ray is on the wrong side of the surface (BRDF is 0 as ray is pointing away from the hemisphere )*/
                            if (surface_bsdf_pdf == 0)
                            {
                                break;
                            }

                            cur_n = hrec.normal;
                            throughput *= bsdf * cos_wi;
                            result += throughput;
                        }
                    }
                    else
                    {
                        // Hit nothing, which means it'll hit envmap. Calculate SH coefficients here
                        const double cosine = std::max(dot(cur_n, unit_vector(r.direction())), 0.0);
                        const double cos_wi = std::max(dot(n, unit_vector(v_direction)), 0.0);
                        for (int l = 0; l < n_bands; ++l)
                        {
                            for (int m = -l; m <= l; ++m)
                            {
                                double phi = std::atan2(r.direction().x(), -r.direction().z());
                                double theta = std::acos(std::clamp(r.direction().y(), -1.0, 1.0));
                                phi = (phi < 0) ? (phi + M_PI * 2) : phi;
                                theta = (theta < 0) ? (theta + M_PI) : theta;
                                const double SH_basis_sample = EstimateSH(l, m, theta, phi);
                                const Vector3f value = (cosine * bsdf * result + v_bsdf * cos_wi) * SH_basis_sample;
                                tri->coeffs[idx][l * (l + 1) + m] += (value);
                            }
                        }
                        break;
                    }
                }
            }
            // Divide the result by number of samples
            const double factor = 4.0 * M_PI / (double)n_samples;
            for (int i = 0; i < n_coeffs; ++i)
                tri->coeffs[idx][i] *= factor;
        }
    });
    tf.dispatch().get();
}

// Here, n_coeffs = n_bands*n_bands and n_samples = sqrt_n_samples*sqrt_n_samples
void path_prt::SH_project_environment()
{
    Li_coeffs = {};
    image* img = dynamic_cast<image_texture*>(scene->env_map->env_map_tex.get())->img.get();
    image_texture* tex = dynamic_cast<image_texture*>(scene->env_map->env_map_tex.get());
    for (int y = 0; y < img->ny; ++y)
    {
        for (int x = 0; x < img->nx; ++x)
        {
            const double theta = M_PI * (y + 0.5) / (double)img->ny;
            const double phi = 2 * M_PI * (x + 0.5) / (double)img->nx;
            Vector3f pixel = tex->value(x, y);
            for (int l = 0; l < n_bands; ++l)
            {
                for (int m = -l; m <= l; ++m)
                {
                    const double SH_basis_sample = EstimateSH(l, m, theta, phi);
                    Li_coeffs[l * (l + 1) + m] += pixel * SH_basis_sample;
                }
            }
        }
    }
    // Divide the result by weight and number of samples
    const double factor = 2.0 * M_PI * M_PI / (img->nx * img->ny);
    for (int i = 0; i < n_coeffs; ++i)
    {
        Li_coeffs[i] *= factor;
    }
}

// TODO
// Convert from recursive to iterative
Vector3f path_prt::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const double &prev_bsdf_pdf, sampler &random_sampler)
{
    hit_record hrec;
    auto &world = scene->world;
    Vector3f result(0.0, 0.0, 0.0);

    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec(hrec);
        if (hrec.mat_ptr->scatter(r, hrec, srec, random_sampler.get3d()))
        {
            auto &tri_coeffs = dynamic_cast<triangle*>(hrec.obj)->coeffs;
            Vector2f &uv = hrec.uv;
            SHCoefficients tr_coeffs = {};

            for (int i = 0; i < n_coeffs; ++i)
            {
                tr_coeffs[i] += (1 - uv.x - uv.y) * tri_coeffs[0][i] + uv.x * tri_coeffs[1][i] + uv.y * tri_coeffs[2][i];
                result += Li_coeffs[i] * tr_coeffs[i];
            }
            return result;
        }
    }
    return scene->env_map->eval(r, hrec, depth);
}

