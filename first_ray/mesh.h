#ifndef MESH_H_
#define MESH_H_

#include "hitable.h"
#include <vector>

class triangle : public hitable
{
public:
    triangle() {}
    triangle(vec3 _v0, vec3 _v1, vec3 _v2, material *mat) : v0(_v0), v1(_v1), v2(_v2), mat_ptr(mat) {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        const float EPSILON = 0.0000001;
        float a, f, u, v;
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec3 h = cross(r.d, edge2);
        a = dot(edge1, h);
        // Check if this ray is parallel to this triangle's plane.
        if (a > -EPSILON && a < EPSILON)
            return false;    
        f = 1.0 / a;
        vec3 s = r.o - v0;
        u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        vec3 q = cross(s, edge1);
        v = f * dot(r.d, q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        float t = f * dot(edge2, q);

        // Check if ray intersected successfully
        if (t > EPSILON)
        {
            rec.t = t;
            rec.p = r.point_at_parameter(t);
            rec.normal = cross(edge1, edge2);
            rec.normal.make_unit_vector();
            rec.u = u;
            rec.v = v;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return false;
    }

    vec3 v0;
    vec3 v1;
    vec3 v2;

    material *mat_ptr;
};

class mesh : public hitable
{
public:
    std::vector<triangle> triangles;


};


#endif