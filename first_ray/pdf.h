#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"
#include "util.h"
#include <algorithm>

inline Vector3f hemisphere_to_cosine_direction(const Vector2f &sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float r = sqrt(r0);
    const float phi = 2 * (float)M_PI * r1;
    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}

inline Vector3f hemisphere_to_cosine_direction(float &theta, float& phi, const Vector2f &sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float r = sqrt(r0);
    theta = asin(r);
    phi = 2 * (float)M_PI * r1;
    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}


inline Vector3f hemisphere_to_cosine_power_direction(const Vector2f& sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float sin_theta = sqrt(1 - r0);
    const float r = 2 * (float)M_PI *r1;

    return Vector3f(sin_theta*cos(r), sin_theta * sin(r), sqrt(r1));
}

<<<<<<< HEAD
inline Vector3f uniform_sample_sphere(const Vector2f &sample) {
    Vector2f u(sample);
=======
inline Vector3f uniform_sample_sphere() {
    Vector2f u(gen_cano_rand(), gen_cano_rand());
>>>>>>> 789dbd6... Refactor viewer into a separate class
    const float  z = 1 - 2 * u[0];
    const float r = std::sqrt(std::max((float)0, (float)1 - z * z));
    const float phi = 2 * (float)M_PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline Vector3f random_to_sphere(const float &radius, const float &distance_squared, const Vector2f& sample)
{
    const float &r1 = sample[0];
    const float &r2 = sample[1];
    float z = 1 + r2 * (sqrt(1 - radius*radius/distance_squared) - 1);
    float phi = 2*(float)M_PI*r1;
    float x = cos(phi)*sqrt(1-z*z);
    float y = sin(phi)*sqrt(1-z*z);

    return Vector3f(x, y, z);
}

class pdf
{
public:
    virtual float value(const hit_record&, const Vector3f&) const = 0;
    virtual Vector3f generate(const Vector2f &sample) const = 0;
    virtual ~pdf() = 0;
};

class cosine_pdf : public pdf
{
public:
    cosine_pdf(const Vector3f &w) { uvw.build_from_w(w); }
    float value(const hit_record &hrec, const Vector3f &direction) const override
    {
        float cosine = std::max<float>(dot(uvw.w(), unit_vector(direction)), 0.0f);
        
        return cosine / (float) (float)M_PI;
    }

    Vector3f generate(const Vector2f &sample) const override
    {
        return uvw.local(hemisphere_to_cosine_direction(sample));
    }

    onb uvw;
};

class cosine_power_pdf : public pdf
{
public:
    cosine_power_pdf(const ray &r_in, const Vector3f &w, const float &specular_exponent) : specular_exponent(specular_exponent), r_in(r_in) { uvw.build_from_w(w); }
    // TODO: modify this stuff - modified phong BRDF paper
    float value(const hit_record &hrec, const Vector3f &direction) const override
    {
        const float cosine = std::max<float>(0.0f, dot(reflect(unit_vector(r_in.direction()), uvw.w()), direction));

        return ((specular_exponent + 1) / (2 * (float)M_PI)) * pow(cosine, specular_exponent);
    }

    Vector3f generate(const Vector2f &sample) const override
    {
        return uvw.local(hemisphere_to_cosine_power_direction(sample));
    }

    onb uvw;
    const float specular_exponent;
    const ray r_in;
};

class constant_pdf : public pdf
{
public:
    constant_pdf(const float &p) : p(p) { }
    float value(const hit_record &hrec, const Vector3f &direction) const override
    {
        return p;
    }
    
    Vector3f generate(const Vector2f &sample) const override
    {
        throw std::runtime_error("Not implemented!");
    }

    const float p;
};

class mixture_pdf : public pdf
{
public:
    mixture_pdf(std::unique_ptr<pdf> p0, std::unique_ptr<pdf> p1, const float &prob_pdf, const float &sample) : prob_pdf(prob_pdf), rand_var(sample)
    {
        p[0] = std::move(p0);
        p[1] = std::move(p1);
    }

    float value(const hit_record &rec, const Vector3f &direction) const override
    {
        return prob_pdf*p[0]->value(rec, direction) + (1 - prob_pdf)*p[1]->value(rec, direction);
    }

    Vector3f generate(const Vector2f &sample) const override
    {
        if (rand_var < prob_pdf)
            return p[0]->generate(sample);

        return p[1]->generate(sample);
    }

    std::unique_ptr<pdf> p[2];
    const float prob_pdf;
    const float rand_var;
};

#endif
