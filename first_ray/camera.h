#ifndef CAMERA_H_
#define CAMERA_H_

#include "ray.h"

Vector3f random_in_unit_disk()
{
    Vector3f p;
    do
    {
        p = 2.0f * Vector3f(drand48(), drand48(), 0.0f) - Vector3f(1.0f, 1.0f, 0.0f);
    } while (dot(p, p) >= 1.0f);

    return p;
}

struct camera
{
    camera(Vector3f lookfrom, Vector3f lookat, Vector3f vup, float vfov, float aspect, float aperture, float focus_dist)
    {
        lens_radius = aperture / 2;
        float theta = vfov * (float)M_PI / 180.0f;
        float half_height = tan(theta / 2);
        float half_width = aspect * half_height;

        origin = lookfrom;
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        lower_left_corner = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
        horizontal = 2*half_width*focus_dist*u;
        vertical = 2*half_height*focus_dist*v;

    }

    ray get_ray(const float &s, const float &t)
    {
        Vector3f rd = lens_radius * random_in_unit_disk();
        Vector3f offset = u * rd.x() + v * rd.y();
        return ray(origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset);
    }

    Vector3f origin;
    Vector3f lower_left_corner;
    Vector3f horizontal;
    Vector3f vertical;
    Vector3f u, v, w;
    float lens_radius;

};


#endif