#ifndef TEXTURE_H_
#define TEXTURE_H_

#include "geometry.h"
#include "hitable.h"
#include <cmath>

class texture
{
public:
    virtual Vector3f value(const hit_record &rec) const = 0;
};

class constant_texture : public texture
{
public:
    constant_texture() {}
    constant_texture(const Vector3f &c) : color(c) {}
    virtual Vector3f value(const hit_record &rec) const
    {
        return color;
    }

    Vector3f color;
};

class checker_texture : public texture
{
public:
    checker_texture() {}
    checker_texture(texture *t0, texture *t1) : even(t0), odd(t1) {}
    virtual Vector3f value(const hit_record &rec) const
    {
        const Vector3f &p = rec.p;
        const float &u = rec.u;
        const float &v = rec.v;

        float sines = sin(10*p.x()) * sin(10*p.y()) * sin(10*p.z());
        if (sines < 0.0f)
            return odd->value(rec);
        else
            return even->value(rec);
    }

    texture *odd;
    texture *even;
};

#endif