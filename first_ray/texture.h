#ifndef TEXTURE_H_
#define TEXTURE_H_

#include "vec3.h"

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

#endif