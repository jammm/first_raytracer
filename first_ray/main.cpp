
#include <iostream>
#include <cstdlib>
#include "hitable_list.h"
#include "sphere.h"
#include "camera.h"
#include <float.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"

float drand48()
{
    return float(rand()) / float(RAND_MAX);
}

vec3 random_in_unit_sphere()
{
    vec3 p;

    do
    {
        p = 2.0f * vec3(drand48(), drand48(), drand48()) - vec3(1, 1, 1);
    } while (p.squared_length() >= 1.0f);

    return p;
}


vec3 color(const ray &r, hitable *world)
{
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec))
    {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5f * color( ray(rec.p, target - rec.p), world);
    }
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    return (1.0f - t)*vec3(1.0f, 1.0f, 1.0f) + t*vec3(0.5f, 0.7f, 1.0f);
}

int main()
{
    const int nx = 400;
    const int ny = 200;
    const int ns = 100;
    const int comp = 3; //RGB
    unsigned char out_image[nx * ny * comp];
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *list[2];
    camera cam;
    list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), 0.5f);
    list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f);
    hitable *world = new hitable_list(list, 2);

    stbi_flip_vertically_on_write(true);

    for (int j = ny-1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++)
        {
            vec3 col(0.0f, 0.0f, 0.0f);
            for (int s = 0; s < ns; s++)
            {
                float u = float(i + float(rand()) / float(RAND_MAX)) / float(nx);
                float v = float(j + float(rand()) / float(RAND_MAX)) / float(ny);
                ray r = cam.get_ray(u, v);
                //vec3 p = r.point_at_parameter(2.0f);
                col += color(r, world);
            }
            col /= float(ns);
            int ir = int(sqrt(col[0]) * 255.99);
            int ig = int(sqrt(col[1]) * 255.99);
            int ib = int(sqrt(col[2]) * 255.99);
            int index = (j * nx + i) * comp;
            out_image[index]     = unsigned char(ir);
            out_image[index + 1] = unsigned char(ig);
            out_image[index + 2] = unsigned char(ib);
            //std::cout << ir << " " << ig << " " << ib << "\n";
        }
        std::cout << ".";
        //std::cout << (float(ny - 1 - j) / float(ny)) * 100.0f << "%\r\r\r\r";
    }

    stbi_write_bmp("out.bmp", nx, ny, comp, (void *)out_image);
    stbi_write_png("out.png", nx, ny, comp, (void *)out_image, nx * comp);
    stbi_write_jpg("out.jpg", nx, ny, comp, (void *)out_image, 100);

    return 0;
}