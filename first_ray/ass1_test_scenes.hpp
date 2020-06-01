#ifndef TEST_SCENES_H
#define TEST_SCENES_H

#include "hitable.h"

#ifdef USE_SSE
#include "triangle_sse.hpp"
#else
#include "triangle.h"
#endif
#include "parallel_bvh.h"
#include "material.h"
#include "camera.h"
#include "point_light.h"
#include "Scene.h"

Scene *bunny_scene(const float &aspect)
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

    static std::unique_ptr<material> lambert = std::make_unique<lambertian>(new constant_texture(Vector3f(1, 1, 1)));

    // Table of objects
    obj objects[] =
    {
        {"data/bunny.obj", lambert.get(), identity},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);

    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    std::vector<hitable *> lights;
    lights.push_back(new point_light(Vector3f(10, 20, 10), new point_light_mat(new constant_texture(Vector3f(1000, 1000, 1000)))));

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

    // Triangle floor
    Vector3f *vertices = new Vector3f[3];
    vertices[0] = Vector3f(-100, 0, -100);
    vertices[1] = Vector3f(0, 0, 100);
    vertices[2] = Vector3f(100, 0, -100);

    Vector3f *normals = new Vector3f[3];
    normals[0] = Vector3f(0, 1, 0);
    normals[1] = Vector3f(0, 1, 0);
    normals[2] = Vector3f(0, 1, 0);

    std::vector<int> indices{ 0, 1, 2 };

    Vector2f *uv = new Vector2f[3];

    auto mesh = std::make_shared<triangle_mesh>(1, 3,
        vertices, indices.data(), normals, uv, std::move(lambert), "floor", true);

    list[i++] = new triangle(std::move(mesh), 0);

    Vector3f lookfrom(0, 5, 15);
    Vector3f lookat(0, 0, 0);
    constexpr float dist_to_focus = 100;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 45.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.2f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *bunny_20_scene(const float &aspect)
{
    hitable **list = new hitable*[2500000];
    int i = 0;

    struct obj
    {
        const char *obj_file;
        material *bsdf;
        Matrix4x4 toWorld;
    };

    static std::unique_ptr<material> lambert = std::make_unique<lambertian>(new constant_texture(Vector3f(1, 1, 1)));

    // Define bunny stuff
    std::vector<Matrix4x4> bunnies(20);
    for (auto &bunny : bunnies)
        bunny.set_identity();

    for (auto it = bunnies.begin() + 10; it != bunnies.cend(); ++it)
        *it = (*it).rotate(110, 0, 1, 0).scale(.6, 1, 1.1);

    // Table of objects
    obj objects[] =
    {
        {"data/bunny.obj"  , lambert.get(), bunnies[0].scale(0.3, 2.0, 0.7).translate(-1, .4, .3).rotate(25, .3, .1, .6)},
        {"data/bunny.obj"  , lambert.get(), bunnies[1].scale(.6, 1.2, .9).translate(7.6, .8, .6)},
        {"data/bunny.obj"  , lambert.get(), bunnies[2].translate(.7, 0, -2).rotate(120, 0, .6, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[3].translate(3.6, 3, -1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[4].translate(-2.4, 2, 3).scale(1, .8, 2)},
        {"data/bunny.obj"  , lambert.get(), bunnies[5].translate(5.5, -.5, 1).scale(1, 2, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[6].rotate(15, 0, 0, 1).translate(-4, -.5, -6).scale(1, 2, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[7].rotate(60, 0, 1, 0).translate(5, .1, 3)},
        {"data/bunny.obj"  , lambert.get(), bunnies[8].translate(-3, .4, 6).rotate(-30, 0, 1, 0)},
        {"data/bunny.obj"  , lambert.get(), bunnies[9].translate(3, 0.5, -2).rotate(180, 0, 1, 0).scale(1.5, 1.5, 1.5)},
        {"data/bunny.obj"  , lambert.get(), bunnies[10].scale(0.3, 2.0, 0.7).translate(-1, .4, .3).rotate(25, .3, .1, .6)},
        {"data/bunny.obj"  , lambert.get(), bunnies[11].scale(.6, 1.2, .9).translate(7.6, .8, .6)},
        {"data/bunny.obj"  , lambert.get(), bunnies[12].translate(.7, 0, -2).rotate(120, 0, .6, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[13].translate(3.6, 3, -1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[14].translate(-2.4, 2, 3).scale(1, .8, 2)},
        {"data/bunny.obj"  , lambert.get(), bunnies[15].translate(5.5, -.5, 1).scale(1, 2, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[16].rotate(15, 0, 0, 1).translate(-4, -.5, -6).scale(1, 2, 1)},
        {"data/bunny.obj"  , lambert.get(), bunnies[17].rotate(60, 0, 1, 0).translate(5, .1, 3)},
        {"data/bunny.obj"  , lambert.get(), bunnies[18].translate(-3, .4, 6).rotate(-30, 0, 1, 0)},
        {"data/bunny.obj"  , lambert.get(), bunnies[19].translate(3, 0.5, -2).rotate(180, 0, 1, 0).scale(1.5, 1.5, 1.5)},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);

    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    std::vector<hitable *> lights;
    lights.push_back(new point_light(Vector3f(10, 20, 10), new point_light_mat(new constant_texture(Vector3f(1000, 1000, 1000)))));

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

    // Triangle floor
    Vector3f *vertices = new Vector3f[3];
    vertices[0] = Vector3f(-100, 0, -100);
    vertices[1] = Vector3f(0, 0, 100);
    vertices[2] = Vector3f(100, 0, -100);

    Vector3f *normals = new Vector3f[3];
    normals[0] = Vector3f(0, 1, 0);
    normals[1] = Vector3f(0, 1, 0);
    normals[2] = Vector3f(0, 1, 0);

    std::vector<int> indices{ 0, 1, 2 };

    Vector2f *uv = new Vector2f[3];

    auto mesh = std::make_shared<triangle_mesh>(1, 3,
        vertices, indices.data(), normals, uv, std::move(lambert), "floor", true);

    list[i++] = new triangle(std::move(mesh), 0);

    Vector3f lookfrom(0, 5, 15);
    Vector3f lookat(0, 0, 0);
    constexpr float dist_to_focus = 100;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 45.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.2f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *teapot_scene(const float &aspect)
{
    hitable **list = new hitable*[10000];
    int i = 0;

    struct obj
    {
        const char *obj_file;
        material *bsdf;
        Matrix4x4 toWorld;
    };

    Matrix4x4 identity;
    identity.set_identity();

    static std::unique_ptr<material> lambert = std::make_unique<lambertian>(new constant_texture(Vector3f(1, 1, 1)));

    // Table of objects
    obj objects[] =
    {
        {"data/teapot.obj", lambert.get(), identity},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);

    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    static std::vector<hitable *> lights;
    lights.push_back(new point_light(Vector3f(10, 10, 10), new point_light_mat(new constant_texture(Vector3f(700, 700, 700)))));

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

    // Triangle floor
    Vector3f *vertices = new Vector3f[3];
    vertices[0] = Vector3f(-10, 0, -10);
    vertices[1] = Vector3f(0, 0, 10);
    vertices[2] = Vector3f(10, 0, -10);

    Vector3f *normals = new Vector3f[3];
    normals[0] = Vector3f(0, 1, 0);
    normals[1] = Vector3f(0, 1, 0);
    normals[2] = Vector3f(0, 1, 0);

    std::vector<int> indices{ 0, 1, 2 };

    Vector2f *uv = new Vector2f[3];

    auto mesh = std::make_shared<triangle_mesh>(1, 3,
        vertices, indices.data(), normals, uv, std::move(lambert), "floor", true);

    list[i++] = new triangle(std::move(mesh), 0);

    Vector3f lookfrom(0, 3, 6);
    Vector3f lookat(0, 0, 0);
    constexpr float dist_to_focus = 100;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 45.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0, 0, 0.2f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

Scene *sponza_scene(const float &aspect)
{
    hitable **list = new hitable*[3000000];
    int i = 0;

    struct obj
    {
        const char *obj_file;
        material *bsdf;
        Matrix4x4 toWorld;
    };

    Matrix4x4 identity;
    identity.set_identity();

    material *lambert = new lambertian(new constant_texture(Vector3f(1, 1, 1)));

    // Table of objects
    obj objects[] =
    {
        {"data/sponza.obj", lambert, identity},
    };
    constexpr int num_objs = sizeof(objects) / sizeof(obj);

    static std::vector<std::vector<std::shared_ptr<hitable>>> triangle_soup;
    std::vector<hitable *> lights;
    lights.push_back(new point_light(Vector3f(0, 10.0, 0), new point_light_mat(new constant_texture(Vector3f(200, 200, 200)))));

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

    Vector3f lookfrom(8, 1.5f, 1);
    Vector3f lookat(0, 2.5f, -1);
    constexpr float dist_to_focus = 100;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 55.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map(std::make_unique<constant_texture>(Vector3f(0.0f, 0.0f, 0.2f))),
        cam, lights
    );
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new hitable_list(std::vector<hitable*>(list, list + i), i);
}

#endif