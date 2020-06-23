#ifndef AARECT_H
#define AARECT_H

#include "ray.h"
#include "hitable.h"

class xy_rect : public hitable
{
public:
	xy_rect() = default;
    xy_rect(const double &x0, const double &x1, const double &y0, const double &y1, const double &k, material *mat) : x0(x0), x1(x1), y0(y0), y1(y1), k(k), mat(mat) {}


    virtual bool bounding_box(double t0, double t1, aabb &box) const
    {
        box = aabb(Vector3f(x0, y0, k - 0.0001), Vector3f(x1, y1, k + 0.0001));
        return true;
    }
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const
    {
        const double t = (k - r.origin().z()) / r.direction().z();
        if (t < t_min || t > t_max)
            return false;
        const double x = r.origin().x() + t*r.direction().x();
        const double y = r.origin().y() + t*r.direction().y();

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


    double x0;
    double x1;
    double y0;
    double y1;
    // z coord. of this XY plane
    double k;
    material *mat;
};

class xz_rect : public hitable
{
public:
    xz_rect() {}
    xz_rect(const double &x0, const double &x1, const double &z0, const double &z1, const double &k, material *mat) : x0(x0), x1(x1), z0(z0), z1(z1), k(k), mat(mat) {}


    bool bounding_box(double t0, double t1, aabb &box) const override
    {
        box = aabb(Vector3f(x0, k - 0.0001, z0), Vector3f(x1, k + 0.0001, z1));
        return true;
    }

    virtual double pdf_direct_sampling(const Vector3f &o, const Vector3f &v) const
    {
        hit_record rec;
        if (this->hit(ray(o, v), 0.001, FLT_MAX, rec))
        {
            double area = (x1 - x0)*(z1 - z0);
            double distance_squared = rec.t*rec.t*v.squared_length();
            double cosine = fabs(dot(v, rec.normal) / v.length());

            return distance_squared / (cosine * area);
        }
        return 0;
    }

    virtual Vector3f random(const Vector3f &o, const Vector2f& sample) const
    {
        Vector3f random_point = Vector3f(x0+sample[0]*(x1-x0), k, z0+sample[1]*(z1-z0));
        return random_point - o;
    }

    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override
    {
        const double t = (k - r.origin().y()) / r.direction().y();
        if (t < t_min || t > t_max)
            return false;
        const double x = r.origin().x() + t*r.direction().x();
        const double z = r.origin().z() + t*r.direction().z();

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


    double x0;
    double x1;
    double z0;
    double z1;
    // y coord. of this XZ plane
    double k;
    material *mat;
};

class yz_rect : public hitable
{
public:
    yz_rect() {}
    yz_rect(const double &y0, const double &y1, const double &z0, const double &z1, const double &k, material *mat) : y0(y0), y1(y1), z0(z0), z1(z1), k(k), mat(mat) {}


    virtual bool bounding_box(double t_min, double t_max, aabb &box) const
    {
        box = aabb(Vector3f(k - 0.0001, y0, z0), Vector3f(k + 0.0001, y1, z1));
        return true;
    }
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const
    {
        const double t = (k - r.origin().x()) / r.direction().x();
        if (t < t_min || t > t_max)
            return false;
        const double y = r.origin().y() + t*r.direction().y();
        const double z = r.origin().z() + t*r.direction().z();

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

    double y0;
    double y1;
    double z0;
    double z1;
    // x coord. of this XY plane
    double k;
    material *mat;
};

#endif
