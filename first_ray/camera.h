#ifndef CAMERA_H_
#define CAMERA_H_

#include "ray.h"
#include "util.h"

struct camera
{
    camera() {}
    camera(Vector3f lookfrom, Vector3f lookat, Vector3f vup, double vfov, double aspect, double aperture, double focus_dist)
    {
        lens_radius = aperture / 2;
        double theta = vfov * (double)M_PI / 180.0;
        double half_height = tan(theta / 2);
        double half_width = aspect * half_height;

        origin = lookfrom;
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        lower_left_corner = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
        horizontal = 2*half_width*focus_dist*u;
        vertical = 2*half_height*focus_dist*v;

        dist = 512.0 / (2 * half_height);

    }

    ray get_ray(const double s, const double t, const Vector2f &sample)
    {
        Vector2f rd = (lens_radius == 0) ? Vector2f(lens_radius, lens_radius) : (lens_radius * random_in_unit_disk(sample));
        Vector3f offset = u * rd[0] + v * rd[1];
        return ray(origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset);
    }

    ray get_ray_as_pinhole(const double s, const double t)
    {
        return ray(origin, lower_left_corner + s * horizontal + t * vertical - origin);
    }

    Vector3f origin;
    Vector3f lower_left_corner;
    Vector3f horizontal;
    Vector3f vertical;
    Vector3f u, v, w;
    double lens_radius;
    double dist;

};


#endif
