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
    checker_texture(texture *t0, texture *t1, const float u_scale_, const float v_scale_) : tex0(t0), tex1(t1), u_scale(u_scale_), v_scale(v_scale_) {}
    virtual Vector3f value(const hit_record &hrec) const
    {
        int x = 2 * modulo((int)(hrec.u * u_scale * 2), 2) - 1,
            y = 2 * modulo((int)(hrec.v * v_scale * 2), 2) - 1;

        if (x * y == 1)
            return tex0->value(hrec);
        else
            return tex1->value(hrec);
    }

    texture *tex0;
    texture *tex1;
    float u_scale, v_scale;
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
        const int nx = img->nx;
        const int ny = img->ny;
        int i = (rec.u)*nx;
        int j = (1 - rec.v)*ny;
		//std::cout << "i,j: " << i << " " << j << std::endl;
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i > nx - 1) i = nx - 1;
        if (j > ny - 1) j = nx - 1;
        float r, g, b;
        const int idx = (j * nx + i) * 3;
        if (img->type == formats::STBI_HDR)
        {
            r = img->dataf[idx];
            g = img->dataf[idx + 1];
            b = img->dataf[idx + 2];
        }
        else
        {
            r = int(img->data[idx]) / 255.0f;
            g = int(img->data[idx + 1]) / 255.0f;
            b = int(img->data[idx + 2]) / 255.0f;
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
