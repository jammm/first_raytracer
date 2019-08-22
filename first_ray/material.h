#ifndef MATERIAL_H
#define MATERIAL_H

#include "hitable.h"
#include "util.h"
#include "texture.h"
#include "onb.h"
#include "pdf.h"
#include "image.h"

inline Vector3f random_in_unit_sphere()
{
    Vector3f p;

    do
    {
        p = 2.0f * Vector3f(gen_cano_rand(), gen_cano_rand(), gen_cano_rand()) - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0f);

    return p;
}

inline Vector3f random_on_unit_sphere()
{
    Vector3f p;

    do
    {
        p = 2.0f * Vector3f(gen_cano_rand(), gen_cano_rand(), gen_cano_rand()) - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0f);

    return unit_vector(p);
}

struct scatter_record
{
    ray specular_ray;
    bool is_specular;
    Vector3f attenuation;
    std::unique_ptr<pdf> pdf_ptr;
};

class material
{
public:
    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const { return false; }
    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const { return Vector3f(0, 0, 0); }
    virtual Vector3f get_albedo(const hit_record& rec) const { return Vector3f(0, 0, 0); }
    virtual Vector3f emitted(const ray &r_in, const hit_record &rec) const { return Vector3f(0, 0, 0); }
};

class lambertian : public material
{
public:
    lambertian(texture *a) : albedo(a) {}

    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const
    {
        srec.is_specular = false;
        srec.attenuation = albedo->value(hrec);
        srec.pdf_ptr = std::make_unique<cosine_pdf>(hrec.normal);
        return true;
    }

    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const
    {
        return albedo->value(rec) / M_PI;
    }

    virtual Vector3f get_albedo(const hit_record& rec) const
    {
        return albedo->value(rec);
    }

    texture *albedo;
};

class modified_phong : public material
{
public:
    modified_phong(texture *a, texture *b, const float &specular_exponent) : diffuse_reflectance(a), specular_reflectance(b), specular_exponent(specular_exponent) {}

    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const
    {
        //const float rand_var = gen_cano_rand();
        //constexpr float specular_chance = 0.5f;

        srec.is_specular = true;

        srec.pdf_ptr = std::make_unique<cosine_power_pdf>(r_in, hrec.normal, specular_exponent);
        srec.attenuation = diffuse_reflectance->value(hrec);
        srec.specular_ray = ray(hrec.p, srec.pdf_ptr->generate());

        return true;
    }

    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const
    {
        const Vector3f reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        const float cos_alpha = dot(reflected, r_in.direction());
        return (diffuse_reflectance->value(rec) / M_PI) + specular_reflectance->value(rec)
                * ((specular_exponent + 2) / (2 * M_PI)) * pow(dot(reflected, wo), specular_exponent);
    }

    texture *diffuse_reflectance;
    texture *specular_reflectance;
    const float specular_exponent;
    Vector3f wo;
};

class metal : public material
{
public:
    metal(const Vector3f &a, float f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1.0f; }
    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const
    {
        Vector3f reflected = reflect(unit_vector(r_in.direction()), hrec.normal);
        srec.specular_ray = ray(hrec.p, reflected + fuzz*random_in_unit_sphere());
        srec.attenuation = albedo;
        srec.is_specular = true;
        srec.pdf_ptr = std::make_unique<constant_pdf>(1.0f);
        return true;
    }
    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const
    {
        return albedo;
    }

    Vector3f albedo;
    float fuzz;
};

static float schlick(const float &cosine, const float &ref_idx)
{
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0)*pow((1 - cosine), 5);
}

static bool refract(const Vector3f &v, const Vector3f &n, float ni_over_nt, Vector3f &refracted)
{
    Vector3f uv = unit_vector(v);
    float dt = dot(uv, n);
    float discriminant = 1.0f - ni_over_nt*ni_over_nt*(1 - dt*dt);
    if (discriminant > 0)
    {
        refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
        return true;
    }
    return false;
}

class dielectric : public material
{
public:
    dielectric(float ri) : ref_idx(ri) {}

    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const
    {
        Vector3f outward_normal;
        Vector3f reflected = reflect(r_in.direction(), hrec.normal);
        float ni_over_nt;
        srec.is_specular = true;
        srec.pdf_ptr = 0;
        srec.attenuation = Vector3f(1.0f, 1.0f, 1.0f);
        Vector3f refracted;
        float reflect_prob;
        float cosine;
        if (dot(r_in.direction(), hrec.normal) > 0)
        {
            outward_normal = -hrec.normal;
            ni_over_nt = ref_idx;
            cosine = ref_idx * dot(r_in.direction(), hrec.normal) / r_in.direction().length();
        }
        else
        {
            outward_normal = hrec.normal;
            ni_over_nt = 1.0f / ref_idx;
            cosine = -dot(r_in.direction(), hrec.normal) / r_in.direction().length();
        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
        {
            reflect_prob = schlick(cosine, ref_idx);
        }
        else
        {
            reflect_prob = 1.0f;
        }
        if (gen_cano_rand() < reflect_prob)
        {
            srec.specular_ray = ray(hrec.p, reflected);
        }
        else
        {
            srec.specular_ray = ray(hrec.p, refracted);
        }

        return true;
    }
    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const
    {
        return Vector3f(1, 1, 1);
    }

    float ref_idx;
};

class diffuse_light : public material
{
public:
    diffuse_light(texture *a) : emit(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec) const { return false; }
    virtual Vector3f emitted(const ray &r_in, const hit_record &rec) const
    {
        if (dot(rec.normal, r_in.direction()) < 0)
            return emit->value(rec);
        else
            return Vector3f(0, 0, 0);
    }
    texture *emit;
};

class environment_map : public material
{
public:
    environment_map(std::string env_map_filename) : env_map_filename(env_map_filename)
    {
        auto env_map_img = std::make_unique<image>(env_map_filename);
        env_map_tex = std::make_unique<image_texture>(env_map_img);
    }
    environment_map(std::unique_ptr<texture> e)
    {
        env_map_tex = std::move(e);
    }
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const { return false; }
    Vector3f eval(const ray& r_in, hit_record rec, const int &depth) const
    {
        const Vector3f direction = unit_vector(r_in.d);
        float phi = std::atan2(direction.x(), -direction.z());
        float theta = acos(direction.y());

        phi = (phi < 0) ? (phi + M_PI * 2) : phi;
        theta = (theta < 0) ? (theta + M_PI) : theta;

        rec.u = phi / (2.0f * M_PI);
        rec.v = theta / (M_PI);

        if (depth == 0)
            return FromSrgb(env_map_tex->value(rec));

        return env_map_tex->value(rec);
    }

    std::string env_map_filename;
    std::unique_ptr<texture> env_map_tex;
};

#endif
