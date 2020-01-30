#include <iostream>
#include "Scene.h"
#include "sphere.h"
#include "triangle.h"
#include "material.h"
#include "camera.h"
#include "image.h"
#include "texture.h"
#include "util.h"
#include "bvh.h"
#include "parallel_bvh.h"
#include "aarect.h"
#include "box.h"
#include "pdf.h"
#include "viewer.h"

// Include renderers
#include "integrator.h"
#include "path.h"
#include "path_prt.h"
#include "ao.h"
#include "debug_renderer.h"

// Other includes
#include <float.h>
#include <chrono>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/DefaultLogger.hpp>
#include <assimp/LogStream.hpp>

inline Vector3f de_nan(const Vector3f &c)
{
    Vector3f temp = c;
    if (!(temp[0] == temp[0])) temp[0] = 0;
    if (!(temp[1] == temp[1])) temp[1] = 0;
    if (!(temp[2] == temp[2])) temp[2] = 0;

    return temp;
}

hitable *random_scene(camera &cam, const float &aspect, std::vector<hitable *> &lights)
{
    /* n == number of spheres */
    constexpr int n = 1100000;
    hitable **list = new hitable*[n + 1];
    //Large sphere's texture can be checkered
    texture *checker = new checker_texture(new constant_texture(Vector3f(0.2f, 0.3f, 0.1f)), new constant_texture(Vector3f(0.99f, 0.99f, 0.99f)));

    list[0] = new sphere(Vector3f(0, -1000, 0), 1000, new lambertian(checker));
    int i = 1;
    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = gen_cano_rand();
            Vector3f center(a + 0.9f * gen_cano_rand(), 0.2f, b + 0.9f * gen_cano_rand());
            if ((center - Vector3f(4, 0.2f, 0)).length() > 0.9f)
            {
                if (choose_mat < 0.8)
                {
                    //diffuse
                    list[i++] = new sphere(center, 0.2f, new lambertian(new constant_texture(Vector3f(gen_cano_rand() * gen_cano_rand(), gen_cano_rand() * gen_cano_rand(), gen_cano_rand() * gen_cano_rand()))));
                }
                else if (choose_mat < 0.95)
                {
                    //metal
                    list[i++] = new sphere(center, 0.2f,
                        new metal(Vector3f(0.5f*(1 + gen_cano_rand()), 0.5f*(1 + gen_cano_rand()), 0.5f*(1 + gen_cano_rand())), 0.5f*gen_cano_rand()));
                }
                else
                {
                    //dielectric
                    list[i++] = new sphere(center, 0.2f, new dielectric(1.5f));
                }
            }
        }
    }

    Vector3f lookfrom(12, 2, 3);
    Vector3f lookat(0, 0, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.001f;
    constexpr float vfov = 40.0f;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable *>(list, list + i), i);
}

hitable *cornell_box(camera &cam, const float &aspect)
{
    hitable **list = new hitable*[8];
    int i = 0;
    material *red = new lambertian(new constant_texture(Vector3f(0.65f, 0.05f, 0.05f)));
    material *white = new lambertian(new constant_texture(Vector3f(0.73f, 0.73f, 0.73f)));
    material *green = new lambertian(new constant_texture(Vector3f(0.12f, 0.45f, 0.15f)));
    material *light = new diffuse_light(new constant_texture(Vector3f(15, 15, 15)));

    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new flip_normals(new xz_rect(213, 343, 227, 332, 554, light));
    list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, white));
    list[i++] = new xz_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    material *glass = new dielectric(1.5);
    list[i++] = new sphere(Vector3f(190, 90, 190), 90, glass);
    list[i++] = new translate(new rotate_y(
        new box(Vector3f(0, 0, 0), Vector3f(165, 330, 165), white), 15), Vector3f(265, 0, 295));

    Vector3f lookfrom(278, 278, -800);
    Vector3f lookat(278, 278, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;

    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    //return new hitable_list(list, i);
    return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
}

Scene* prt_test(const float& aspect)
{
    hitable** list = new hitable * [20000];
    int i = 0;

    std::vector<hitable*> lights;
    //lights.push_back(new sphere(Vector3f(0, 0, 0), 1.0f, lightt));
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("cube/plane.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    Vector3f lookfrom(25, 25, 200.0f);
    Vector3f lookat(25, 25, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        new hitable_list(std::vector<hitable*>(list, list + i), i),
        new environment_map("data/test.hdr"),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *furnace_test_scene(const float &aspect)
{
    hitable **list = new hitable*[20000];
    int i = 0;

    //Vector3f reflectance(1.0f, 1.0f, 1.0f);
    //material *specular = new modified_phong(new constant_texture(Vector3f(0.0f, 0.0f, 0.0f)),
    //                                    new constant_texture(reflectance), 100.0f);
    //material *lambert = new lambertian(new constant_texture(Vector3f(1, 1, 1)));
    //material *mirror = new metal(Vector3f(1, 1, 1), 0.0f);
    //material *lightt = new diffuse_light(new constant_texture(Vector3f(1, 1, 1)));

    //list[i++] = new sphere(Vector3f(0, 0, 0), 1.0f, lightt);
    //list[i++] = new sphere(Vector3f(0, 0, 0), 0.3f, lambert);
    //list[i++] = new sphere(Vector3f(0.30f, 50, 0), 50.0f, mirror);

    std::vector<hitable*> lights;
    //lights.push_back(new sphere(Vector3f(0, 0, 0), 1.0f, lightt));
    
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("cube/teapot.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }
    
    Vector3f lookfrom(0, 125, 175.0f);
    Vector3f lookat(0, 50, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map("data/ennis.hdr"),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *cornell_box_obj(const float &aspect)
{
    hitable **list = new hitable*[300];
    int i = 0;
    std::vector<hitable*> lights;
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("CornellBox/CornellBox-MIS-Test.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    Vector3f lookfrom(0, 1, 3.9f);
    Vector3f lookat(0, 1, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;
    camera cam;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    //return new hitable_list(std::vector<hitable *>(list, list + i), i);
    //return new bvh_node(list, i, 0.0f, 0.0f);
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.0f))),
        cam, lights
    );
}

Scene* cornell_box_ao(const float& aspect)
{
    hitable** list = new hitable * [300];
    int i = 0;
    std::vector<hitable*> lights;
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("CornellBox/CornellBox-Original.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    Vector3f lookfrom(0, 1, 3.9f);
    Vector3f lookat(0, 1, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;
    camera cam;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(1.0f, 1.0f, 1.0f))),
        cam, lights
    );
}

hitable *veach_mis(camera &cam, const float &aspect, std::vector<hitable *> &lights)
{
    hitable **list = new hitable*[300];
    int i = 0;
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("veach_mi/veach_mi.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    list[i++] = new sphere(Vector3f(10, 10, 4), 0.5f, new diffuse_light(new constant_texture(Vector3f(800, 800, 800))));
    list[i++] = new sphere(Vector3f(-1.25f, 0, 0), 0.1f, new diffuse_light(new constant_texture(Vector3f(100, 100, 100))));
    list[i++] = new sphere(Vector3f(-3.75f, 0, 0), 0.03333f, new diffuse_light(new constant_texture(Vector3f(901.803f, 901.803f, 901.803f))));
    list[i++] = new sphere(Vector3f(1.25f, 0, 0), 0.3f, new diffuse_light(new constant_texture(Vector3f(11.1111f, 11.1111f, 11.1111f))));
    list[i++] = new sphere(Vector3f(3.75f, 0, 0), 0.9f, new diffuse_light(new constant_texture(Vector3f(1.23457, 1.23457, 1.23457))));

    lights.push_back(new sphere(Vector3f(10, 10, 4), 0.5f, new diffuse_light(new constant_texture(Vector3f(800, 800, 800)))));
    lights.push_back(new sphere(Vector3f(-1.25f, 0, 0), 0.1f, new diffuse_light(new constant_texture(Vector3f(100, 100, 100)))));
    lights.push_back(new sphere(Vector3f(-3.75f, 0, 0), 0.03333f, new diffuse_light(new constant_texture(Vector3f(901.803f, 901.803f, 901.803f)))));
    lights.push_back(new sphere(Vector3f(1.25f, 0, 0), 0.3f, new diffuse_light(new constant_texture(Vector3f(11.1111f, 11.1111f, 11.1111f)))));
    lights.push_back(new sphere(Vector3f(3.75f, 0, 0), 0.9f, new diffuse_light(new constant_texture(Vector3f(1.23457, 1.23457, 1.23457)))));

    Vector3f lookfrom(0, 2, 15);
    Vector3f lookat(0, -2, 2.5);
    constexpr float dist_to_focus = 50.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 28.0f;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new hitable_list(std::vector<hitable *>(list, list + i), i);
    //return new bvh_node(list, i, 0.0f, 0.0f);
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
}

int main(int argc, const char **argv)
{
    constexpr int nx = 1024;
    constexpr int ny = 768;
    int ns = 100;
    constexpr int comp = 3; //RGB
    //out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);
    viewer film_viewer(nx, ny, ns, comp);

    /* Parse command line args */
    if (argc == 3)
    {
        const std::string arg = argv[1];
        if (arg == "--ns")
        {
            ns = std::atoi(argv[2]);
        }
    }

    std::cout<<"Resolution: "<<nx<<"x"<<ny<<std::endl;
    std::cout<<"Setting number of samples to "<<ns<<std::endl;

    // Start performance timer
    std::chrono::high_resolution_clock::time_point t11 = std::chrono::high_resolution_clock::now();

    // Initialize scene
    //std::unique_ptr<hitable> world(cornell_box_obj(cam, float(nx) / float(ny), lights));
    std::unique_ptr<Scene> scene(furnace_test_scene(float(nx) / float(ny)));

    std::chrono::high_resolution_clock::time_point t22 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_spann = std::chrono::duration_cast<std::chrono::duration<double>>(t22 - t11);
    std::cout << "\nBVH construction took me " << time_spann.count() << " seconds.";

    // Use the renderer specified in template parameter
    renderer<path> render;

    render.Render(scene.get(), film_viewer);

    // Save film to file(s) (currently JPG,PNG and PFM)
    film_viewer.save_and_destroy();

    return 0;
}
