#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"
#include <algorithm>

class pdf
{
public:
    virtual float value(const Vector3f &direction) const = 0;
    virtual Vector3f generate() const = 0;
};

class cosine_pdf : public pdf
{
public:
    cosine_pdf(const Vector3f &w) { uvw.build_from_w(w); }
    virtual float value(const Vector3f &direction) const
    {
        float cosine = dot(unit_vector(direction), uvw.w());
        
        return std::max<float>(0, cosine / (float) M_PI);
    }
    virtual Vector3f generate() const
    {
        return uvw.local(random_cosine_direction());
    }

    onb uvw;
};

// hitable_pdf is used to generate random directions and to generate a pdf for the corresponding hitable object
class hitable_pdf : public pdf
{
public:
    hitable_pdf(hitable *p, const Vector3f &origin) : ptr(p), o(origin) {}
    virtual float value(const Vector3f &direction) const
    {
        return ptr->pdf_value(o, direction);
    }
    virtual Vector3f generate() const
    {
        return ptr->random(o);
    }

    Vector3f o;
    hitable *ptr;
};

class mixture_pdf : public pdf
{
public:
    mixture_pdf(pdf *p0, pdf *p1)
    {
        p[0] = p0;
        p[1] = p1;
    }

    virtual float value(const Vector3f &direction) const
    {
        return 0.5f*p[0]->value(direction) + 0.5f*p[1]->value(direction);
    }

    virtual Vector3f generate() const
    {
        if (drand48() < 0.5)
            return p[0]->generate();

        return p[1]->generate();
    }

    pdf *p[2];
};

#endif