#ifndef SPHERE_H_
#define SPHERE_H_

#include <memory>
#include "hitable.h"
#include "pdf.h"

class sphere : public hitable
{
public:
    sphere() {}
    sphere(const Vector3f &cen, const double &r, material *mat) : center(cen), radius(r), mat_ptr(mat) {}

    bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override;
    bool bounding_box(double t0, double t1, aabb &b) const override;
    double pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const override;
    Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f& sample) const override;

    Vector3f center;
    double radius;
    material *mat_ptr;
};

bool sphere::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    Vector3f oc = r.origin() - center;
    const double a = dot(r.direction(), r.direction());
    const double b = dot(oc, r.direction());
    const double c = dot(oc, oc) - radius * radius;

    double discriminant = b*b - a*c;
    
    if (discriminant >= 0.0)
    {
        discriminant = sqrt(discriminant);
        double t = (-b - discriminant) / a;
        if (t < t_min) t = (-b + discriminant) / a;

        if (t < t_min || t > t_max)
            return false;

        rec.t = t;
        rec.p = r.point_at_parameter(rec.t);
        rec.normal = (rec.p - center) / radius;
        if ((r.origin() - center).squared_length() < radius * radius)
            rec.normal = -rec.normal;
        rec.mat_ptr = mat_ptr;
        get_sphere_uv((rec.p - center), rec.u, rec.v);
        return true;
    }
    return false;
};

bool sphere::bounding_box(double t0, double t1, aabb &box) const
{
	box = aabb(center - Vector3f(radius, radius, radius), center + Vector3f(radius, radius, radius));
	return true;
}

double sphere::pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const
{
    const Vector3f o = ray(lrec.p, to_light).point_at_parameter(-lrec.t);
    const Vector3f direction = center - o;
    const double distance_squared = direction.squared_length();
    const double radius_squared = radius * radius;
    // If origin is inside sphere, return uniform sphere sampling PDF
    if (distance_squared <= radius_squared)
        return 1 / (4 * M_PI * radius * radius);

    const double cos_theta_max = sqrt(1 - radius_squared / (center-o).squared_length());
    const double solid_angle = 2 * M_PI * (1-cos_theta_max);
    return 1 / solid_angle;
}

Vector3f sphere::sample_direct(hit_record &rec, const Vector3f &o, const Vector2f& sample) const
{
    const Vector3f direction = center - o;
    const double distance_squared = direction.squared_length();
    onb uvw;
    uvw.build_from_w(direction);
    // If origin is inside sphere, do uniform sphere sampling instead
    if (distance_squared <= radius * radius)
    {
        Vector3f p = center + radius * uniform_sample_sphere(sample);
        rec.p = p;
        rec.mat_ptr = mat_ptr;
        rec.normal = unit_vector(center - p);
		rec.obj = (hitable *)this;

        return p - o;
    }
    return uvw.fromLocal(random_to_sphere(radius, distance_squared, sample));
}

#endif
