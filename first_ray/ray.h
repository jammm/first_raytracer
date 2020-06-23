#ifndef RAY_H_
#define RAY_H_

#include "geometry.h"

struct ray
{
    ray() {};
    ray(const Vector3f &a, const Vector3f &b) { o = a; d = b; }
    Vector3f origin() const { return o; }
    Vector3f direction() const { return d; }

    Vector3f point_at_parameter(const double &t) const { return o + d * t; }

    Vector3f o;
    Vector3f d;
};

struct ray4
{
    union
    {
        __m128 ox4, oy4, oz4;
        double o[4];
    };
    union
    {
        __m128 dx4, dy4, dz4;
        double d[4];
    };
    union {
        __m128 t4;
        double t[4];
    };
};

#endif
