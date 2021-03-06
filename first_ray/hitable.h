#ifndef HITABLE_H_
#define HITABLE_H_

#include "ray.h"
#ifdef USE_SSE
#include "aabb_sse.hpp"
#else
#include "aabb.h"
#endif
#include <cfloat>

class material;
class hitable;

inline void get_sphere_uv(const Vector3f &p, double &u, double &v)
{
    double phi = atan2(p.z(), p.x());
    double theta = asin(p.y());
    u = 1 - (phi + M_PI) / (2 * M_PI);
    v = (theta + M_PI / 2) / M_PI;
}

struct alignas(16) hit_record
{
    double t;
	// Global texture coordinates
    double u;
	// Local barycentric coordinates
	Vector2f uv;
    Vector3f p;
    Vector3f wi;
    Vector3f normal;
    material *mat_ptr;
	hitable *obj;
    double v;
};

class hitable
{
public:
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const = 0;
	virtual bool bounding_box(double t0, double t1, aabb &box) const = 0;
    virtual double pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const { return 0.0; }
    virtual Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f &sample) const { return Vector3f(1, 0, 0); }
    virtual ~hitable() = 0;

    // Some global members for SSE stuff
    static const __m128 EPS4;
    static const __m128 MINUSEPS4;
    static const __m128 ONE4;
};

class flip_normals : public hitable
{
public:
    flip_normals(hitable *ptr) : ptr(ptr) {}
    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override
    {
        if (ptr->hit(r, t_min, t_max, rec))
        {
            rec.normal = -rec.normal;
            return true;
        }
        return false;
    }
    bool bounding_box(double t0, double t1, aabb &box) const override
    {
        return ptr->bounding_box(t0, t1, box);
    }

    hitable *ptr;
};

class translate : public hitable
{
public:
    translate(hitable *ptr, const Vector3f &offset) : ptr(ptr), offset(offset) {}
    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override
    {
        ray moved_r(r.origin() - offset, r.direction());
        if (ptr->hit(moved_r, t_min, t_max, rec))
        {
            rec.p += offset;
            return true;
        }
        return false;
    }
    bool bounding_box(double t0, double t1, aabb &box) const override
    {
        if (ptr->bounding_box(t0, t1, box))
        {
            box = aabb(box.min + offset, box.max + offset);
            return true;
        }
        return false;
    }

    hitable *ptr;
    Vector3f offset;
};

class rotate_y : public hitable
{
public:
    rotate_y(hitable *ptr, const double &angle) : ptr(ptr)
    {
        double radians = (M_PI / 180.0) * angle;
        sin_theta = sin(radians);
        cos_theta = cos(radians);
        hasbox = ptr->bounding_box(0, 1, bbox);
        Vector3f min(FLT_MAX, FLT_MAX, FLT_MAX);
        Vector3f max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                {
                    double x = i * bbox.max.x() + (1 - i) * bbox.min.x();
                    double y = j * bbox.max.y() + (1 - j) * bbox.min.y();
                    double z = k * bbox.max.z() + (1 - k) * bbox.min.z();
                    double newx = cos_theta*x + sin_theta*x;
                    double newz = -sin_theta*x + cos_theta*z;
                    Vector3f tester(newx, y, newz);
                    for (int c = 0; c < 3; ++c)
                    {
                        if (tester[c] > max[c])
                            max[c] = tester[c];
                        if (tester[c] < min[c])
                            min[c] = tester[c];
                    }
                }
        bbox = aabb(min, max);
    }

    bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override
    {
        Vector3f origin = r.origin();
        Vector3f direction = r.direction();
        origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
        origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];
        direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
        direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];
        ray rotated_r(origin, direction);
        if (ptr->hit(rotated_r, t_min, t_max, rec))
        {
            Vector3f p = rec.p;
            Vector3f normal = rec.normal;
            p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
            p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];
            normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
            normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];
            rec.p = p;
            rec.normal = normal;
            return true;
        }
        return false;
    }

    bool bounding_box(double t0, double t1, aabb &box) const override
    {
        box = bbox;
        return hasbox;
    }

    hitable *ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

#endif
