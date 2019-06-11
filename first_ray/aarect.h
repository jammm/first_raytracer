#ifndef AARECT_H
#define AARECT_H

#include "ray.h"
#include "hitable.h"

class xy_rect : public hitable
{
public:
    xy_rect() {}
    xy_rect(const float &x0, const float &x1, const float &y0, const float &y1, const float &k, material *mat) : x0(x0), x1(x1), y0(y0), y1(y1), k(k), mat(mat) {}


    virtual bool bounding_box(float t0, float t1, aabb &box) const
    {
        box = aabb(Vector3f(x0, y0, k - 0.0001), Vector3f(x1, y1, k + 0.0001));
        return true;
    }
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        const float t = (k - r.origin().z()) / r.direction().z();
        if (t < t_min || t > t_max)
            return false;
        const float x = r.origin().x() + t*r.direction().x();
        const float y = r.origin().y() + t*r.direction().y();

        // check if ray lies outside the rect
        if (x < x0 || y < y0 || x > x1 || y > y1)
            return false;

        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(0, 0, 1);

        return true;
    }


    float x0;
    float x1;
    float y0;
    float y1;
    // z coord. of this XY plane
    float k;
    material *mat;
};

class xz_rect : public hitable
{
public:
    xz_rect() {}
    xz_rect(const float &x0, const float &x1, const float &z0, const float &z1, const float &k, material *mat) : x0(x0), x1(x1), z0(z0), z1(z1), k(k), mat(mat) {}


    virtual bool bounding_box(float t0, float t1, aabb &box) const
    {
        box = aabb(Vector3f(x0, k - 0.0001, z0), Vector3f(x1, k + 0.0001, z1));
        return true;
    }

    virtual float pdf_direct_sampling(const Vector3f &o, const Vector3f &v) const
    {
        hit_record rec;
        if (this->hit(ray(o, v), 0.001, FLT_MAX, rec))
        {
            float area = (x1 - x0)*(z1 - z0);
            float distance_squared = rec.t*rec.t*v.squared_length();
            float cosine = fabs(dot(v, rec.normal) / v.length());

            return distance_squared / (cosine * area);
        }
        return 0;
    }
    virtual Vector3f random(const Vector3f &o) const
    {
        Vector3f random_point = Vector3f(x0+gen_cano_rand()*(x1-x0), k, z0+gen_cano_rand()*(z1-z0));
        return random_point - o;
    }

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        const float t = (k - r.origin().y()) / r.direction().y();
        if (t < t_min || t > t_max)
            return false;
        const float x = r.origin().x() + t*r.direction().x();
        const float z = r.origin().z() + t*r.direction().z();

        // check if ray lies outside the rect
        if (x < x0 || z < z0 || x > x1 || z > z1)
            return false;

        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(0, 1, 0);

        return true;
    }


    float x0;
    float x1;
    float z0;
    float z1;
    // y coord. of this XZ plane
    float k;
    material *mat;
};

class yz_rect : public hitable
{
public:
    yz_rect() {}
    yz_rect(const float &y0, const float &y1, const float &z0, const float &z1, const float &k, material *mat) : z0(z0), z1(z1), y0(y0), y1(y1), k(k), mat(mat) {}


    virtual bool bounding_box(float t_min, float t_max, aabb &box) const
    {
        box = aabb(Vector3f(k - 0.0001, y0, z0), Vector3f(k + 0.0001, y1, z1));
        return true;
    }
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        const float t = (k - r.origin().x()) / r.direction().x();
        if (t < t_min || t > t_max)
            return false;
        const float y = r.origin().y() + t*r.direction().y();
        const float z = r.origin().z() + t*r.direction().z();

        // check if ray lies outside the rect
        if (z < z0 || y < y0 || z > z1 || y > y1)
            return false;

        rec.u = (y - y0) / (y1 - y0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(1, 0, 0);

        return true;
    }

    float y0;
    float y1;
    float z0;
    float z1;
    // x coord. of this XY plane
    float k;
    material *mat;
};

#endif
