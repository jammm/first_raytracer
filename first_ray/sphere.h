#ifndef SPHERE_H_
#define SPHERE_H_

#include "hitable.h"

class sphere : public hitable
{
public:
    sphere() {}
    sphere(const Vector3f &cen, const float &r, material *mat) : center(cen), radius(r), mat_ptr(mat) {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;

    Vector3f center;
    float radius;
    material *mat_ptr;
};

bool sphere::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    Vector3f oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc, oc) - radius * radius;

    float discriminant = b*b - a*c;
    
    if (discriminant > 0.0f)
    {
        float temp = (-b - sqrt(b*b - a*c)) / a;
        if (temp < t_max && temp > t_min)
        {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            get_sphere_uv((rec.p - center), rec.u, rec.v);
            return true;
        }
        temp = (-b + sqrt(b*b - a*c)) / a;
        if (temp < t_max && temp > t_min)
        {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            get_sphere_uv((rec.p - center), rec.u, rec.v);
            return true;
        }
    }
    return false;
};

#endif