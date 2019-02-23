#ifndef TEXTURE_H_
#define TEXTURE_H_

#include "vec3.h"
#include <cmath>

class texture
{
public:
    virtual vec3 value(const float &u, const float &v, const vec3 &p) const = 0;
};

class constant_texture : public texture
{
public:
    constant_texture() {}
    constant_texture(const vec3 &c) : color(c) {}
    virtual vec3 value(const float &u, const float &v, const vec3 &p) const
    {
        return color;
    }

    vec3 color;
};

class checker_texture : public texture
{
public:
    checker_texture() {}
    checker_texture(texture *t0, texture *t1) : even(t0), odd(t1) {}
    virtual vec3 value(const float &u, const float &v, const vec3 &p) const
    {
        float sines = sin(10*p.x()) * sin(10*p.y()) * sin(10*p.z());
        if (sines < 0.0f)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

    texture *odd;
    texture *even;
};

#endif