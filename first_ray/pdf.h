#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"
#include "util.h"
#include <algorithm>

inline Vector3f hemisphere_to_cosine_direction()
{
    const float r0 = gen_cano_rand(), r1 = gen_cano_rand();
    const float r = sqrt(r0);
    const float phi = 2 * M_PI * r1;
    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}

inline Vector3f hemisphere_to_cosine_power_direction()
{
    const float r0 = gen_cano_rand(), r1 = gen_cano_rand();
    const float sin_theta = sqrt(1 - r0);
    const float r = 2 * M_PI*r1;

    return Vector3f(sin_theta*cos(r), sin_theta * sin(r), sqrt(r1));
}

inline Vector3f random_to_sphere(const float &radius, const float &distance_squared)
{
    float r1 = gen_cano_rand();
    float r2 = gen_cano_rand();
    float z = 1 + r2 * (sqrt(1 - radius*radius/distance_squared) - 1);
    float phi = 2*M_PI*r1;
    float x = cos(phi)*sqrt(1-z*z);
    float y = sin(phi)*sqrt(1-z*z);

    return Vector3f(x, y, z);
}

class pdf
{
public:
    virtual float value(const hit_record &lrec, const Vector3f &to_light) const = 0;
    virtual Vector3f generate() const = 0;
};

class cosine_pdf : public pdf
{
public:
    cosine_pdf(const Vector3f &w) { uvw.build_from_w(w); }
    virtual float value(const hit_record &hrec, const Vector3f &direction) const
    {
        float cosine = std::max<float>(dot(uvw.w(), unit_vector(direction)), 0.0f);
        
        return cosine / (float) M_PI;
    }
    virtual Vector3f generate() const
    {
        return uvw.local(hemisphere_to_cosine_direction());
    }

    onb uvw;
};

class cosine_power_pdf : public pdf
{
public:
    cosine_power_pdf(const ray &r_in, const Vector3f &w, const float &specular_exponent) : r_in(r_in), specular_exponent(specular_exponent) { uvw.build_from_w(w); }
    // TODO: modify this stuff - modified phong BRDF paper
    virtual float value(const hit_record &hrec, const Vector3f &direction) const
    {
        float cosine = std::max<float>(0.0f, dot(reflect(unit_vector(r_in.direction()), uvw.w()), direction));

        return ((specular_exponent + 1) / (2 * M_PI)) * pow(cosine, specular_exponent);
    }
    virtual Vector3f generate() const
    {
        return uvw.local(hemisphere_to_cosine_power_direction());
    }

    onb uvw;
    const float specular_exponent;
    const ray r_in;
};

// hitable_pdf is used to generate random directions and to generate a pdf for the corresponding hitable object
// currently used *only* for sampling on a light source
class hitable_pdf : public pdf
{
public:
    hitable_pdf(hitable *p, const hit_record &hrec) : ptr(p), origin(hrec.p) {}
    virtual float value(const hit_record &lrec, const Vector3f &to_light) const
    {
        return ptr->pdf_direct_sampling(lrec, to_light);
    }
    virtual Vector3f generate(hit_record &rec) const
    {
        return ptr->sample_direct(rec, origin);
    }

    // hitable object on which sample is generated
    hitable *ptr;
    // origin from where direction towards light is sampled
    Vector3f origin;
};

class mixture_pdf : public pdf
{
public:
    mixture_pdf(std::unique_ptr<pdf> p0, std::unique_ptr<pdf> p1, const float &prob_pdf, const float &rand_var = gen_cano_rand()) : prob_pdf(prob_pdf), rand_var(rand_var)
    {
        p[0] = std::move(p0);
        p[1] = std::move(p1);
    }

    virtual float value(const hit_record &rec, const Vector3f &direction) const
    {
        return prob_pdf*p[0]->value(rec, direction) + (1 - prob_pdf)*p[1]->value(rec, direction);
    }

    virtual Vector3f generate() const
    {
        if (rand_var < prob_pdf)
            return p[0]->generate();

        return p[1]->generate();
    }

    std::unique_ptr<pdf> p[2];
    const float prob_pdf;
    const float rand_var;
};

#endif
