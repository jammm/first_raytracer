#ifndef TEXTURE_H_
#define TEXTURE_H_

#include <cmath>
#include "geometry.h"
#include "hitable.h"
#include "image.h"
#include "util.h"

class texture
{
public:
    virtual Vector3f value(const hit_record &rec) const = 0;
    virtual ~texture() = 0;
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
    checker_texture(texture *t0, texture *t1) : odd(t1), even(t0) {}
    virtual Vector3f value(const hit_record &hrec) const
    {
        int x = 2 * modulo((int)(hrec.uv[0] * 2), 2) - 1,
            y = 2 * modulo((int)(hrec.uv[1] * 2), 2) - 1;

        if (x * y == 1)
            return even->value(hrec);
        else
            return odd->value(hrec);
    }

    texture *odd;
    texture *even;
};

class image_texture : public texture
{
public:
    image_texture() {}
    image_texture(std::unique_ptr<image> img_)
    {
        img = std::move(img_);
    }

    Vector3f value(const hit_record &rec) const override
    {
        const int &nx = img->nx;
        const int &ny = img->ny;
        int i = (rec.u)*nx;
        int j = (rec.v)*ny;
		//std::cout << "i,j: " << i << " " << j << std::endl;
        if (i < 0) return Vector3f(0, 0, 0);
        if (j < 0) return Vector3f(0, 0, 0);
        if (i > nx - 1) return Vector3f(0, 0, 0);
        if (j > ny - 1) return Vector3f(0, 0, 0);
        float r, g, b;
        if (img->type == formats::STBI_HDR)
        {
            r = img->dataf[3 * i + 3 * nx * j];
            g = img->dataf[3 * i + 3 * nx * j + 1];
            b = img->dataf[3 * i + 3 * nx * j + 2];
        }
        else
        {
            r = int(img->data[3 * i + 3 * nx * j]) / 255.0f;
            g = int(img->data[3 * i + 3 * nx * j + 1]) / 255.0f;
            b = int(img->data[3 * i + 3 * nx * j + 2]) / 255.0f;
        }

        return Vector3f(r, g, b);
    }

    Vector3f value(const int& x, const int& y) const
    {
        const int& nx = img->nx;
        float r, g, b;
        if (img->type == formats::STBI_HDR)
        {
            r = img->dataf[3 * x + 3 * nx * y];
            g = img->dataf[3 * x + 3 * nx * y + 1];
            b = img->dataf[3 * x + 3 * nx * y + 2];
        }
        else
        {
            r = int(img->data[3 * x + 3 * nx * y]) / 255.0f;
            g = int(img->data[3 * x + 3 * nx * y + 1]) / 255.0f;
            b = int(img->data[3 * x + 3 * nx * y + 2]) / 255.0f;
        }

        return Vector3f(r, g, b);
    }

    std::unique_ptr<image> img;
};

#endif
