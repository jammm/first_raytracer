#ifndef RAY_H_
#define RAY_H_

#include "vec3.h"

struct ray
{
    ray() {};
    ray(const vec3 &a, const vec3 &b) { o = a; d = b; }
    vec3 origin() const { return o; }
    vec3 direction() const { return d; }

    vec3 point_at_parameter(const float &t) const { return o + d * t; }

    vec3 o;
    vec3 d;
};

#endif