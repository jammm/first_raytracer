#ifndef MATERIAL_H
#define MATERIAL_H

#include "hitable.h"
#include "util.h"
#include "texture.h"

Vector3f random_in_unit_sphere()
{
    Vector3f p;

    do
    {
        p = 2.0f * Vector3f(drand48(), drand48(), drand48()) - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0f);

    return p;
}

class material
{
public:
    virtual bool scatter(const ray &r_in, const hit_record &rec, Vector3f &attenuation, ray &scattered) const = 0;
    virtual Vector3f emitted(const hit_record &rec) const { return Vector3f(0, 0, 0); }
};

class lambertian : public material
{
public:
    lambertian(texture *a) : albedo(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, Vector3f &attenuation, ray &scattered) const
    {
        Vector3f target = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p);
        hit_record record = rec;
        record.u = 0;
        record.v = 0;
        attenuation = albedo->value(rec);
        return true;
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
    virtual bool scatter(const ray &r_in, const hit_record &rec, Vector3f &attenuation, ray &scattered) const
    {
        Vector3f reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
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

    virtual bool scatter(const ray &r_in, const hit_record &rec, Vector3f &attenuation, ray &scattered) const
    {
        Vector3f outward_normal;
        Vector3f reflected = reflect(r_in.direction(), rec.normal);
        float ni_over_nt;
        attenuation = Vector3f(1.0f, 1.0f, 1.0f);
        Vector3f refracted;
        float reflect_prob;
        float cosine;
        if (dot(r_in.direction(), rec.normal) > 0)
        {
            outward_normal = -rec.normal;
            ni_over_nt = ref_idx;
            cosine = ref_idx * dot(r_in.direction(), rec.normal) / r_in.direction().length();
        }
        else
        {
            outward_normal = rec.normal;
            ni_over_nt = 1.0f / ref_idx;
            cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
        {
            reflect_prob = schlick(cosine, ref_idx);
        }
        else
        {
            scattered = ray(rec.p, reflected);
            reflect_prob = 1.0f;
        }
        if (drand48() < reflect_prob)
        {
            scattered = ray(rec.p, reflected);
        }
        else
        {
            scattered = ray(rec.p, refracted);
        }

        return true;
    }


    float ref_idx;
};

class diffuse_light : public material
{
public:
    diffuse_light(texture *a) : emit(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, Vector3f &attenuation, ray &scattered) const { return false; }
    virtual Vector3f emitted(const hit_record &rec) const
    { 
        return emit->value(rec);
    }
    texture *emit;
};

#endif