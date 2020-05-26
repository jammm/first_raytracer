#include <iostream>
#include "Scene.h"
#include "sphere.h"
#include "triangle.h"
#include "point_light.h"
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
#include "sampler.h"

// Test scenes
#include "ass1_test_scenes.hpp"

// Include renderers
#include "integrator.h"
#include "path.h"
#include "path_prt.h"
#include "ao.h"
#include "debug_renderer.h"
#include "pssmlt.h"
#include "assignment_1.hpp"

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
    texture *checker = new checker_texture(new constant_texture(Vector3f(0.2f, 0.3f, 0.1f)), new constant_texture(Vector3f(0.99f, 0.99f, 0.99f)), 1, 1);
    sampler s;

    list[0] = new sphere(Vector3f(0, -1000, 0), 1000, new lambertian(checker));
    int i = 1;
    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = s.get1d();
            Vector3f center(a + 0.9f * s.get1d(), 0.2f, b + 0.9f * s.get1d());
            if ((center - Vector3f(4, 0.2f, 0)).length() > 0.9f)
            {
                if (choose_mat < 0.8)
                {
                    //diffuse
                    list[i++] = new sphere(center, 0.2f, new lambertian(new constant_texture(Vector3f(s.get1d() * s.get1d(), s.get1d() * s.get1d(), s.get1d() * s.get1d()))));
                }
                else if (choose_mat < 0.95)
                {
                    //metal
                    list[i++] = new sphere(center, 0.2f,
                        new metal(Vector3f(0.5f*(1 + s.get1d()), 0.5f*(1 + s.get1d()), 0.5f*(1 + s.get1d())), 0.5f*s.get1d()));
                }
                else
                {
                    //dielectric
                    list[i++] = new sphere(center, 0.2f, new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))));
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
    material *glass = new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)));
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
    //material *specular = new modified_phong(new constant_texture(Vector3f(0, 0, 0)),
    //                                    new constant_texture(reflectance), 100.0f);
    //material *lambert = new lambertian(new constant_texture(Vector3f(1, 1, 1)));
    //material *mirror = new metal(Vector3f(1, 1, 1), 0.0f);
    //material *lightt = new diffuse_light(new constant_texture(Vector3f(1, 1, 1)));

    //list[i++] = new sphere(Vector3f(0, 0, 0), 1.0f, lightt);
    //list[i++] = new sphere(Vector3f(0, 0, 0), 0.3f, lambert);
    //list[i++] = new sphere(Vector3f(0.30f, 50, 0), 50.0f, mirror);

    std::vector<hitable*> lights;
    //lights.push_back(new sphere(Vector3f(0, 0, 0), 1.0f, lightt));
#define TEST_MATERIAL new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))
    Matrix4x4 I;
    I.set_identity();
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("cube/teapot.obj", I, TEST_MATERIAL, lights);


    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }
    
    Vector3f lookfrom(0, 125, 175.0f);
    Vector3f lookat(0, 50, 0);
    constexpr float dist_to_focus = 300.0f;
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
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("CornellBox/CornellBox-Original.obj", lights);
    lights.push_back(new point_light(Vector3f(0, 1.5f, 1.9f), new point_light_mat(new constant_texture(Vector3f(100, 100, 100)))));

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

Scene *veach_door_scene(const float &aspect)
{
    hitable **list = new hitable*[500000];
    int i = 0;

    struct obj
    {
        const char *obj_file;
        material *bsdf;
        Matrix4x4 toWorld;
    };

    Matrix4x4 identity;
    identity.set_identity();
    
    Matrix4x4 pot2_Mesh000_toWorld(-0.0757886f, 0, -0.0468591f, -1.95645f, 0, 0.0891049f, 0, 0.648205f, 0.0468591, 0, -0.0757886f, -1.77687f, 0, 0, 0, 1);
    Matrix4x4 pot2_Mesh009_toWorld(-0.0757886f, 0, -0.0468591f, -1.95645f, 0, 0.0891049f, 0, 0.648205f, 0.0468591, 0, -0.0757886f, -1.77687f, 0, 0, 0, 1);
    Matrix4x4 pot3_Mesh000_toWorld(-0.0891049f, 0, 7.7898e-009f, -1.95645f, 0, 0.0891049f, 0, 0.648205f, -7.7898e-009f, 0, -0.0891049f, -2.67687f, 0, 0, 0, 1);
    Matrix4x4 pot3_Mesh009_toWorld(-0.0891049f, 0, 7.7898e-009f, -1.95645f, 0, 0.0891049f, 0, 0.648205f, -7.7898e-009f, 0, -0.0891049f, -2.67687f, 0, 0, 0, 1);
    Matrix4x4 dielectric_pot_Mesh009_toWorld(-0.0837611f, 0, 0.0303939f, -1.95645f, 0, 0.0891049f, 0, 0.651268f, -0.0303939f, 0, -0.0837611f, -3.57687f, 0, 0, 0, 1);
    Matrix4x4 dielectric_pot_Mesh000_toWorld(-0.083761f, 0, 0.0303938f, -1.95645f, 0, 0.0891049f, 0, 0.651268f, -0.0303938f, 0, -0.083761f, -3.57687f, 0, 0, 0, 1);
    Matrix4x4 Mesh011_toWorld(1.8f, 0, 0, 2.3f, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    Matrix4x4 Mesh014_toWorld(1.8f, 0, 0, 2.3f, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);

#define POT2_MATERIAL new rough_conductor(0.15f, 1.0f, Vector3f(1.65746, 0.880369, 0.521229), Vector3f(9.22387, 6.26952, 4.837), new constant_texture(Vector3f(1, 1, 1)), "ggx")
#define FLOOR_MATERIAL new rough_conductor(0.1f, 1.0f, Vector3f(1.65746, 0.880369, 0.521229), Vector3f(9.22387, 6.26952, 4.837), \
    new checker_texture(new constant_texture(Vector3f(0.8, 0.8, 0.8)), new constant_texture(Vector3f(0.2, 0.2, 0.2)), 20, 80), "ggx")
#define DOORHANDLE_MATERIAL new rough_conductor(0.25f, 1.0f, Vector3f(1.65746, 0.880369, 0.521229), Vector3f(9.22387, 6.26952, 4.837), new constant_texture(Vector3f(1, 1, 1)), "beckmann")
#define HINGE_MATERIAL new rough_conductor(0.1f, 1.0f, Vector3f(1.65746, 0.880369, 0.521229), Vector3f(9.22387, 6.26952, 4.837), new constant_texture(Vector3f(1, 1, 1)), "ggx")

    // Table of objects
    obj objects[] =
    {
        {"veach_ajar/models/light.obj"  , nullptr                                                                                                                  , identity},
        {"veach_ajar/models/Mesh000.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , pot3_Mesh000_toWorld},
        {"veach_ajar/models/Mesh009.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , pot3_Mesh009_toWorld},
        {"veach_ajar/models/Mesh000.obj", POT2_MATERIAL                                                                                                            , pot2_Mesh000_toWorld},
        {"veach_ajar/models/Mesh009.obj", POT2_MATERIAL                                                                                                            , pot2_Mesh009_toWorld},                                                                                  
        {"veach_ajar/models/Mesh000.obj", new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , dielectric_pot_Mesh000_toWorld},
        {"veach_ajar/models/Mesh009.obj", new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , dielectric_pot_Mesh009_toWorld},
        {"veach_ajar/models/Mesh001.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , identity},
        {"veach_ajar/models/Mesh002.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , identity},
        {"veach_ajar/models/Mesh003.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , identity},
        {"veach_ajar/models/Mesh004.obj", new lambertian(new image_texture(std::make_unique<image>("veach_ajar/textures/Good Textures_005844.jpg")))               , identity},
        {"veach_ajar/models/Mesh005.obj", new lambertian(new constant_texture(Vector3f(0.247059f, 0.168627f, 0.0901961f)))                                         , identity},
        {"veach_ajar/models/Mesh006.obj", new lambertian(new image_texture(std::make_unique<image>("veach_ajar/textures/cherry-wood-texture.jpg")))                , identity},
        {"veach_ajar/models/Mesh007.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , identity},
        {"veach_ajar/models/Mesh008.obj", new lambertian(new image_texture(std::make_unique<image>("veach_ajar/textures/landscape-with-a-lake.jpg")))              , identity},
        {"veach_ajar/models/Mesh010.obj", HINGE_MATERIAL                                                                                                           , identity},
        {"veach_ajar/models/Mesh011.obj", FLOOR_MATERIAL                                                                                                           , Mesh011_toWorld},
        {"veach_ajar/models/Mesh012.obj", HINGE_MATERIAL                                                                                                           , identity},
        {"veach_ajar/models/Mesh013.obj", new lambertian(new constant_texture(Vector3f(0.258824f, 0.207843f, 0.145098f)))                                          , identity},
        {"veach_ajar/models/Mesh014.obj", new lambertian(new constant_texture(Vector3f(0.8f, 0.8f, 0.8f)))                                                         , Mesh014_toWorld},
        {"veach_ajar/models/Mesh015.obj", DOORHANDLE_MATERIAL                                                                                                      , identity},
        {"veach_ajar/models/Mesh016.obj", HINGE_MATERIAL                                                                                                           , identity},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);
    
    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    std::vector<hitable*> lights;

    std::cout << std::endl;
    for (int j = 0; j < num_objs; ++j)
    {
        std::cout << "\nProcessing mesh #" << j;
        triangle_soup.push_back(create_triangle_mesh(objects[j].obj_file, objects[j].toWorld, objects[j].bsdf, lights));
    }

    for (auto &mesh : triangle_soup)
    {
        for (auto triangle : mesh)
        {
            list[i++] = triangle.get();
        }
    }

    Vector3f lookfrom(4.05402f, 1.61647f, -2.30652f);
    Vector3f lookat(3.064f, 1.58417f, -2.44373f);
    constexpr float dist_to_focus = 100;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 36.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.0f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *glass_of_water(const float &aspect)
{
    hitable **list = new hitable*[1000000];
    int i = 0;

    struct obj
    {
        const char *obj_file;
        material *bsdf;
        Matrix4x4 toWorld;
    };

    Matrix4x4 identity;
    identity.set_identity();

#define BACKDROP_MATERIAL new rough_conductor(0.1f, 1.0f, Vector3f(4.27751f, 3.51315f, 2.76113f), Vector3f(3.49118f, 2.88936f, 3.1117f), new constant_texture(Vector3f(0.578596f, 0.578596f, 0.578596f)), "ggx")
#define FLOOR_MAT new rough_conductor(0.1f, 1.0f, Vector3f(1.65746f, 0.880369f, 0.521229f), Vector3f(9.22387f, 6.26952f, 4.837f), new constant_texture(Vector3f(0.578596f, 0.578596f, 0.578596f)), "ggx")

    // Table of objects
    obj objects[] =
    {
        {"glass-of-water/models/light.obj"  , nullptr                                                                                                                    , identity},
        {"glass-of-water/models/Mesh008.obj", new dielectric(1.33f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh005.obj", new dielectric(1.33f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh004.obj", new dielectric(1.31f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh003.obj", new dielectric(0.763f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , identity},
        {"glass-of-water/models/Mesh010.obj", new dielectric(1.31f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh009.obj", new dielectric(0.763f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , identity},
        {"glass-of-water/models/Mesh012.obj", new dielectric(0.763f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , identity},
        {"glass-of-water/models/Mesh014.obj", new dielectric(1.31f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh015.obj", new dielectric(0.763f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , identity},
        {"glass-of-water/models/Mesh006.obj", new dielectric(1.31f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh002.obj", new dielectric(0.763f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f))) , identity},
        {"glass-of-water/models/Mesh001.obj", new dielectric(1.31f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
        {"glass-of-water/models/Mesh007.obj", BACKDROP_MATERIAL                                                                                                          , identity},
        {"glass-of-water/models/Mesh011.obj", FLOOR_MAT                                                                                                                  , identity},
        {"glass-of-water/models/Mesh013.obj", new dielectric(1.5f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))   , identity},
        {"glass-of-water/models/Mesh000.obj", new dielectric(1.33f, new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)), new constant_texture(Vector3f(1.0f, 1.0f, 1.0f)))  , identity},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);
    
    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    std::vector<hitable*> lights;

    std::cout << std::endl;
    for (int j = 0; j < num_objs; ++j)
    {
        std::cout << "\nProcessing mesh #" << j;
        triangle_soup.push_back(create_triangle_mesh(objects[j].obj_file, objects[j].toWorld, objects[j].bsdf, lights));
    }

    for (auto &mesh : triangle_soup)
    {
        for (auto triangle : mesh)
        {
            list[i++] = triangle.get();
        }
    }

    Vector3f lookfrom(-0.0893585f, 2.69412f, 25.6726f);
    Vector3f lookat(-0.085623f, 2.70998f, 24.6727f);
    constexpr float dist_to_focus = 300.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 21.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.0f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

int main(int argc, const char **argv)
{
    constexpr int nx = 512;
    constexpr int ny = 512;
    int ns = 1;
    constexpr int comp = 3; //RGB
    //out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);

    /* Parse command line args */
    if (argc == 3)
    {
        const std::string arg = argv[1];
        if (arg == "--ns")
        {
            ns = std::atoi(argv[2]);
        }
    }
	viewer film_viewer(nx, ny, ns, comp);
    
	std::cout<<"Resolution: "<<nx<<"x"<<ny<<std::endl;
    std::cout<<"Setting number of samples to "<<ns<<std::endl;

    // Start performance timer
    std::chrono::high_resolution_clock::time_point t11 = std::chrono::high_resolution_clock::now();

    // Initialize scene
    //std::unique_ptr<hitable> world(cornell_box_obj(cam, float(nx) / float(ny), lights));
    std::unique_ptr<Scene> scene(bunny_20_scene(float(nx) / float(ny)));

    std::chrono::high_resolution_clock::time_point t22 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_spann = std::chrono::duration_cast<std::chrono::duration<double>>(t22 - t11);
    std::cout << "\nBVH construction took me " << time_spann.count() << " seconds.";

    // Use the renderer specified in template parameter
    renderer<ass_1_renderer> render;

    render.Render(scene.get(), film_viewer);

    // Save film to file(s) (currently JPG,PNG and PFM)
    film_viewer.save_and_destroy();

    return 0;
}
