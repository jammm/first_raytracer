#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"

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
        if (cosine > 0)
            return cosine / M_PI;
        else
            return 0;
    }
    virtual Vector3f generate() const
    {
        return uvw.local(random_cosine_direction());
    }

    onb uvw;
};

#endif