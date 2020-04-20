#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H

#include "hitable.h"
#include "material.h"

class point_light : public hitable
{
public:
    point_light(const Vector3f p_light_, point_light_mat *mat_) : p_light(p_light_), mat(mat_)
    {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &hrec) const
    {
        // It's point light, so we will never intersect it
        return false;
    }

    virtual bool bounding_box(float t0, float t1, aabb &b) const
    {
        // Not used
        return false;
    }

    virtual float pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const
    {
        return 1.0f;
    }
    virtual Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f& sample) const
    {
        rec.t = 1.0f;
        rec.p = p_light;
        rec.normal = unit_vector(o - p_light);
        rec.mat_ptr = mat;
        rec.u = 0;
        rec.v = 0;
        rec.uv.x = 0;
        rec.uv.y = 0;
		rec.obj = (hitable *)this;
        rec.wi = rec.normal;

        return p_light - o;
    }

    // Position of light
    Vector3f p_light;
    point_light_mat* mat;
};

#endif