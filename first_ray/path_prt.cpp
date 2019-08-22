#include "path_prt.h"
#include "prt.h"
#include "triangle.h"

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
                float theta = samples[i].theta;
                float phi = samples[i].phi;
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
    for (auto& i : world->list)
    {
        triangle* tri = dynamic_cast<triangle*>(i);
        if (tri == nullptr) continue;
        for (int idx = 0; idx < 3; ++idx)
        {
            const Vector3f& v = tri->mesh->vertices[tri->V[idx]];
            const Vector3f& n = tri->mesh->normals[tri->V[idx]];
            const Point2f uv = tri->mesh->uv[tri->V[idx]];
            for (int i = 0; i < n_samples; ++i)
            {
                const Vector3f& direction = samples[i].direction;
                const ray r(v, direction);
                hit_record rec;
                if (!scene->world->hit(r, EPSILON, FLT_MAX, rec))
                {
                    const float cosine = std::max(dot(n, direction), 0.0f);
                    if (cosine == 0.0f) continue;
                    rec.p = v;
                    rec.u = uv.x;
                    rec.v = uv.y;
                    const Vector3f albedo = tri->mat_ptr->get_albedo(rec);
                    float theta = samples[i].theta;
                    float phi = samples[i].phi;
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
            }
        }
    }
}

// Here, n_coeffs = n_bands*n_bands and n_samples = sqrt_n_samples*sqrt_n_samples
void path_prt::SH_project_environment()
{
    Li_coeffs.reset(new Vector3f[n_coeffs]);
    std::memset(Li_coeffs.get(), 0, n_coeffs * sizeof(Vector3f));

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
            }

            for (int i = 0; i < n_coeffs; ++i)
            {
                result += Li_coeffs[i] * tr_coeffs[i];
            }
            return result;
        }
    }
    return scene->env_map->eval(r, hrec, depth);
}

