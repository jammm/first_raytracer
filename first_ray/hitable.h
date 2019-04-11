#ifndef HITABLE_H_
#define HITABLE_H_

#include "ray.h"
#include "aabb.h"

class material;

inline void get_sphere_uv(const Vector3f &p, float &u, float &v)
{
    float phi = atan2(p.z(), p.x());
    float theta = asin(p.y());
    u = 1 - (phi + M_PI) / (2 * M_PI);
    v = (theta + M_PI / 2) / M_PI;
}

struct hit_record
{
    float t;
    float u;
    float v;
    Vector3f p;
    Vector3f normal;
    material *mat_ptr;
};

class hitable
{
public:
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const = 0;
	virtual bool bounding_box(float t0, float t1, aabb &box) const = 0;
};

#endif