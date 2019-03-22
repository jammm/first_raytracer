#ifndef RAY_H_
#define RAY_H_

#include "geometry.h"

struct ray
{
    ray() {};
    ray(const Vector3f &a, const Vector3f &b) { o = a; d = b; }
    Vector3f origin() const { return o; }
    Vector3f direction() const { return d; }

    Vector3f point_at_parameter(const float &t) const { return o + d * t; }

    Vector3f o;
    Vector3f d;
};

#endif