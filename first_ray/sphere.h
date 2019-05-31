#ifndef SPHERE_H_
#define SPHERE_H_

#include "hitable.h"
#include "pdf.h"

class sphere : public hitable
{
public:
    sphere() {}
    sphere(const Vector3f &cen, const float &r, material *mat) : center(cen), radius(r), mat_ptr(mat) {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;
    virtual bool bounding_box(float t0, float t1, aabb &b) const;
    virtual float pdf_value(const Vector3f &o, const Vector3f &v) const;
    virtual Vector3f random(const Vector3f &o) const;

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

bool sphere::bounding_box(float t0, float t1, aabb &box) const
{
	box = aabb(center - Vector3f(radius, radius, radius), center + Vector3f(radius, radius, radius));
	return true;
}

float sphere::pdf_value(const Vector3f &o, const Vector3f &v) const
{
    hit_record rec;
    if (this->hit(ray(o, v), 0.001f, FLT_MAX, rec))
    {
        float cos_theta_max = sqrt(1 - radius*radius / (center-o).squared_length());
        float solid_angle = 2 * M_PI * (1-cos_theta_max);
        return 1 / solid_angle;
    }
    else
        return 0;
}

Vector3f sphere::random(const Vector3f &o) const
{
    const Vector3f direction = center - o;
    const float distance_squared = direction.squared_length();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}

#endif
