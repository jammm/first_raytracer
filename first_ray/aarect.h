#ifndef AARECT_H
#define AARECT_H

#include "ray.h"
#include "hitable.h"

class xyrect : public hitable
{
public:
    xyrect() {}
    xyrect(const float &x0, const float &x1, const float &y0, const float &y1, material *mat) : x0(x0), x1(x1), y0(y0), y1(y1), mat(mat) {}


    virtual bool bounding_box(float t0, float t1, aabb &box)
    {
        box = aabb(Vector3f(x0, y0, k - 0.0001), Vector3f(x1, y1, k + 0.0001));
        return true;
    }
    virtual bool hit(const ray &r, float t0, float t1, hit_record &rec) const
    {
        const float t = (k - r.origin().z()) / r.direction().z();
        const float x = r.origin().x() + t * r.direction().x;
        const float y = r.origin().y() + t * r.direction().y;

        // check if ray lies outside the rect
        if (x < x0 || y < y0 || x > x1 || y > y1)
            return false;

        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(0, 0, 1);
    }


    float x0;
    float x1;
    float y0;
    float y1;
    material *mat;

    // z coord. of this XY plane
    float k;
};

class xzrect : public hitable
{
public:
    xzrect() {}
    xzrect(const float &x0, const float &x1, const float &z0, const float &z1, material *mat) : x0(x0), x1(x1), z0(z0), z1(z1), mat(mat) {}


    virtual bool bounding_box(float t0, float t1, aabb &box)
    {
        box = aabb(Vector3f(x0, z0, k - 0.0001), Vector3f(x1, z1, k + 0.0001));
        return true;
    }
    virtual bool hit(const ray &r, float t0, float t1, hit_record &rec) const
    {
        const float t = (k - r.origin().z()) / r.direction().z();
        const float x = r.origin().x() + t * r.direction().x;
        const float z = r.origin().z() + t * r.direction().z;

        // check if raz lies outside the rect
        if (x < x0 || z < z0 || x > x1 || z > z1)
            return false;

        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(0, 0, 1);
    }


    float x0;
    float x1;
    float z0;
    float z1;
    material *mat;

    // z coord. of this Xz plane
    float k;
};

class yzrect : public hitable
{
public:
    yzrect() {}
    yzrect(const float &y0, const float &y1, const float &z0, const float &z1, material *mat) : z0(z0), z1(z1), y0(y0), y1(y1), mat(mat) {}


    virtual bool bounding_boz(float t0, float t1, aabb &boz)
    {
        boz = aabb(Vector3f(z0, y0, k - 0.0001), Vector3f(z1, y1, k + 0.0001));
        return true;
    }
    virtual bool hit(const ray &r, float t0, float t1, hit_record &rec) const
    {
        const float t = (k - r.origin().z()) / r.direction().z();
        const float z = r.origin().z() + t * r.direction().z;
        const float y = r.origin().y() + t * r.direction().y;

        // check if ray lies outside the rect
        if (z < z0 || y < y0 || z > z1 || y > y1)
            return false;

        rec.u = (z - z0) / (z1 - z0);
        rec.v = (y - y0) / (y1 - y0);
        rec.t = t;
        rec.mat_ptr = mat;
        rec.p = r.point_at_parameter(t);
        rec.normal = Vector3f(0, 0, 1);
    }

    float y0;
    float y1;
    float z0;
    float z1;
    material *mat;

    // z coord. of this zY plane
    float k;
};

#endif