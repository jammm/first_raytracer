
#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include "material.h"
#include "camera.h"
#include "util.h"
#include <float.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"




hitable *random_scene()
{
    int n = 500;
    hitable **list = new hitable*[n + 1];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
    int i = 1;
    for (int a = -11;a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = drand48();
            vec3 center(a + 0.9* drand48(), 0.2, b + 0.9 * drand48());
            if ((center - vec3(4, 0.2, 0)).length() > 0.9)
            {
                if (choose_mat < 0.8)
                {
                    //diffuse
                    list[i++] = new sphere(center, 0.2, new lambertian(vec3(drand48() * drand48(), drand48() * drand48(), drand48() * drand48())));
                }
                else if (choose_mat < 0.95)
                {
                    //metal
                    list[i++] = new sphere(center, 0.2,
                        new metal(vec3(0.5* (1 + drand48()), 0.5 * (1 + drand48()), 0.5 * (1 + drand48())), 0.5*drand48()));
                }
                else
                {
                    //dielectric
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hitable_list(list, i);
}

vec3 color(const ray &r, hitable *world, int depth)
{
    hit_record rec;
    if (world->hit(r, 0.001f, FLT_MAX, rec))
    {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        {
            return attenuation * color(scattered, world, depth + 1);
        }
        else
        {
            return vec3(0.0f, 0.0f, 0.0f);
        }
    }
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    return (1.0f - t)*vec3(1.0f, 1.0f, 1.0f) + t*vec3(0.5f, 0.7f, 1.0f);
}

int main()
{
    const int nx = 1200;
    const int ny = 800;
    const int ns = 100;
    const int comp = 3; //RGB
    unsigned char *out_image = new unsigned char[nx * ny * comp];
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *list[5];
    vec3 lookfrom(13, 2, 3);
    vec3 lookat(0, 0, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.1f;
    camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus);
    float R = cos(M_PI / 4);

    //list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), 0.5f, new lambertian(vec3(0.1f, 0.2f, 0.5f)));
    //list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f, new lambertian(vec3(0.8f, 0.8f, 0.0f)));
    //list[2] = new sphere(vec3(1.0f, 0, -1.0f), 0.5f, new metal(vec3(0.8f, 0.6f, 0.2f), 0.5f));
    //list[3] = new sphere(vec3(-1.0f, 0.0f, -1.0f), 0.5f, new dielectric(1.5f));
    //list[4] = new sphere(vec3(-1.0f, 0.0f, -1.0f), -0.45f, new dielectric(1.5f));
    //hitable *world = new hitable_list(list, 5);
    hitable *world = random_scene();

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
                col += color(r, world, 0);
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

    delete out_image;

    return 0;
}