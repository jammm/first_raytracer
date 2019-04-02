#ifndef TEXTURE_H_
#define TEXTURE_H_

#include "geometry.h"
#include "hitable.h"
#include "image.h"
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

class image_texture : public texture
{
public:
    image_texture() {}
    image_texture(std::unique_ptr<image> &img) : img(move(img)) {}

    virtual Vector3f value(const hit_record &rec) const
    {
        int comp = img->nn;
        const int &nx = img->nx;
        const int &ny = img->ny;
        int i = (rec.u)*nx;
        int j = (1 - rec.v)*ny - 0.001;
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i > nx - 1) i = nx - 1;
        if (j > ny - 1) j = ny - 1;
        float r = int(img->data[3 * i + 3 * nx*j]) / 255.0;
        float g = int(img->data[3 * i + 3 * nx*j + 1]) / 255.0;
        float b = int(img->data[3 * i + 3 * nx*j + 2]) / 255.0;

        return Vector3f(r, g, b);
    }

    std::unique_ptr<image> img;
};

#endif