
#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include <float.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"


vec3 color(const ray &r, hitable *world)
{
    hit_record rec;
    if (world->hit(r, 0.0, FLT_MAX, rec))
    {
        return 0.5f * vec3(rec.normal.x() + 1.0f, rec.normal.y() + 1.0f, rec.normal.z() + 1.0f);
    }
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    return (1.0f - t)*vec3(1.0f, 1.0f, 1.0f) + t*vec3(0.5f, 0.7f, 1.0f);
}

int main()
{
    const int nx = 400;
    const int ny = 200;
    const int comp = 3; //RGB
    unsigned char out_image[nx * ny * comp];
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    vec3 lower_left_corner(-2.0f, -1.0f, -1.0f);
    vec3 horizontal(4.0f, 0.0f, 0.0f);
    vec3 vertical(0.0f, 2.0f, 0.0f);
    vec3 origin(0.0f, 0.0f, 0.0f);
    hitable *list[2];
    list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), 0.5f);
    list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f);
    hitable *world = new hitable_list(list, 2);

    stbi_flip_vertically_on_write(true);

    for (int j = ny-1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++)
        {
            float u = float(i) / float(nx);
            float v = float(j) / float(ny);
            ray r(origin, lower_left_corner + u*horizontal + v*vertical);
            vec3 p = r.point_at_parameter(2.0f);
            vec3 col = color(r, world);

            int ir = int(col[0] * 255.99);
            int ig = int(col[1] * 255.99);
            int ib = int(col[2] * 255.99);
            int index = (j * nx + i) * comp;
            out_image[index]     = unsigned char(ir);
            out_image[index + 1] = unsigned char(ig);
            out_image[index + 2] = unsigned char(ib);
            //std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }

    stbi_write_bmp("out.bmp", nx, ny, comp, (void *)out_image);
    stbi_write_png("out.png", nx, ny, comp, (void *)out_image, nx * comp);
    stbi_write_jpg("out.jpg", nx, ny, comp, (void *)out_image, 100);

    return 0;
}