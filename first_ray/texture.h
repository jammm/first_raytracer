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
    checker_texture(texture *t0, texture *t1, const double u_scale_, const double v_scale_) : tex0(t0), tex1(t1), u_scale(u_scale_), v_scale(v_scale_) {}
    virtual Vector3f value(const hit_record &hrec) const
    {
        int x = 2 * modulo((int)(hrec.u * u_scale * 2), 2) - 1,
            y = 2 * modulo((int)(hrec.v * v_scale * 2), 2) - 1;

        if (x * y == 1)
            return tex1->value(hrec);
        else
            return tex0->value(hrec);
    }

    texture *tex0;
    texture *tex1;
    double u_scale, v_scale;
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
        int j = (rec.v)*ny;
		//std::cout << "i,j: " << i << " " << j << std::endl;
        if (i < 0 || i > nx)
            i = modulo(i, nx);
        if (j < 0 || j > ny)
            j = modulo(j, ny);
        if (i == nx) i = nx - 1;
        if (j == ny) j = ny - 1;
        double r, g, b;
        const int idx = (j * nx + i) * 3;
        if (img->type == formats::STBI_HDR)
        {
            r = img->dataf[idx];
            g = img->dataf[idx + 1];
            b = img->dataf[idx + 2];
        }
        else
        {
            r = FromSrgb(img->data[idx] / 255.0);
            g = FromSrgb(img->data[idx + 1] / 255.0);
            b = FromSrgb(img->data[idx + 2] / 255.0);
        }

        return Vector3f(r, g, b);
    }

    Vector3f value(const int& x, const int& y) const
    {
        const int& nx = img->nx;
        double r, g, b;
        if (img->type == formats::STBI_HDR)
        {
            r = img->dataf[3 * x + 3 * nx * y];
            g = img->dataf[3 * x + 3 * nx * y + 1];
            b = img->dataf[3 * x + 3 * nx * y + 2];
        }
        else
        {
            r = int(img->data[3 * x + 3 * nx * y]) / 255.0;
            g = int(img->data[3 * x + 3 * nx * y + 1]) / 255.0;
            b = int(img->data[3 * x + 3 * nx * y + 2]) / 255.0;
        }

        return Vector3f(r, g, b);
    }

    std::unique_ptr<image> img;
};

#endif
