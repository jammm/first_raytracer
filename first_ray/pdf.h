#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"
#include "util.h"
#include <algorithm>

inline Vector3f random_cosine_direction()
{
    float r0 = gen_cano_rand(), r1 = gen_cano_rand();
    float r = sqrt(r0);
    float phi = 2 * M_PI * r1;
    float x = r * cos(phi);
    float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
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
        return uvw.local(random_cosine_direction());
    }

    onb uvw;
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
    mixture_pdf(pdf *p0, pdf *p1)
    {
        p[0] = p0;
        p[1] = p1;
    }

    virtual float value(const hit_record &rec, const Vector3f &direction) const
    {
        return 0.5f*p[0]->value(rec, direction) + 0.5f*p[1]->value(rec, direction);
    }

    virtual Vector3f generate() const
    {
        if (gen_cano_rand() < 0.5)
            return p[0]->generate();

        return p[1]->generate();
    }

    pdf *p[2];
};

#endif
