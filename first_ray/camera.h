#ifndef CAMERA_H_
#define CAMERA_H_

#include "ray.h"

struct camera
{
    camera()
    {
        lower_left_corner = vec3(-2.0f, -1.0f, -1.0f);
        horizontal = vec3(4.0f, 0.0f, 0.0f);
        vertical = vec3(0.0f, 2.0f, 0.0f);
        origin = vec3(0.0f, 0.0f, 0.0f);
    }

    ray get_ray(const float &u, const float &v) { return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin); }

    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;

};


#endif