#ifndef MATERIAL_H
#define MATERIAL_H

#include "hitable.h"
#include "util.h"
#include "texture.h"
#include "onb.h"
#include "pdf.h"

Vector3f random_in_unit_sphere()
{
    Vector3f p;

    do
    {
        p = 2.0f * Vector3f(drand48(), drand48(), drand48()) - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0f);

    return p;
}

Vector3f random_on_unit_sphere()
{
    Vector3f p;

    do
    {
        p = 2.0f * Vector3f(drand48(), drand48(), drand48()) - Vector3f(1, 1, 1);
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
    virtual float scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const { return false; }
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

    virtual float scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const
    {
        float cosine = dot(rec.normal, unit_vector(scattered.direction()));
        if (cosine < 0.0f)
            cosine = 0.0f;

        return cosine / M_PI;
    }

    texture *albedo;
};

Vector3f reflect(const Vector3f &v, const Vector3f &n)
{
    return v - 2 * dot(v, n)*n;
}

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
        srec.pdf_ptr = nullptr;
        return true;
    }

    Vector3f albedo;
    float fuzz;
};

float schlick(const float &cosine, const float &ref_idx)
{
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0)*pow((1 - cosine), 5);
}

bool refract(const Vector3f &v, const Vector3f &n, float ni_over_nt, Vector3f &refracted)
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
        if (drand48() < reflect_prob)
        {
            srec.specular_ray = ray(hrec.p, reflected);
        }
        else
        {
            srec.specular_ray = ray(hrec.p, refracted);
        }

        return true;
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

#endif