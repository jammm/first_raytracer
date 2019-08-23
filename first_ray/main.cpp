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
#include "path.h"
#include "path_prt.h"
#include "debug_renderer.h"
#include <float.h>
#include <taskflow/taskflow.hpp>
#include <chrono>

#if defined(_WIN32) || defined(__linux__)
#define GLEW_STATIC
#include <GL/glew.h>
#define _CRT_SECURE_NO_DEPRECATE
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#endif
#include <GLFW/glfw3.h>

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

void glfw_error_callback(int, const char* err_str)
{
    std::cout << "GLFW Error: " << err_str << std::endl;
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

    Vector3f reflectance(1.0f, 1.0f, 1.0f);
    material* specular = new modified_phong(new constant_texture(Vector3f(0.0f, 0.0f, 0.0f)),
        new constant_texture(reflectance), 100.0f);
    material* lambert = new lambertian(new constant_texture(Vector3f(1, 1, 1)));
    material* mirror = new metal(Vector3f(1, 1, 1), 0.0f);
    material* lightt = new diffuse_light(new constant_texture(Vector3f(1, 1, 1)));

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

    Vector3f reflectance(1.0f, 1.0f, 1.0f);
    material *specular = new modified_phong(new constant_texture(Vector3f(0.0f, 0.0f, 0.0f)),
                                        new constant_texture(reflectance), 100.0f);
    material *lambert = new lambertian(new constant_texture(Vector3f(1, 1, 1)));
    material *mirror = new metal(Vector3f(1, 1, 1), 0.0f);
    material *lightt = new diffuse_light(new constant_texture(Vector3f(1, 1, 1)));

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
    

    Vector3f lookfrom(50, 125, 175.0f);
    Vector3f lookat(0, 50, 0);
    constexpr float dist_to_focus = 10.0f;
    constexpr float aperture = 0.0f;
    constexpr float vfov = 40.0f;
    camera cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new Scene(
        parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f),
        new environment_map("data/small_empty_house_4k.hdr"),
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

void background_thread(const std::shared_future<void> &future, GLubyte *out_image, GLFWwindow* window, int nx, int ny, bool &to_exit)
{
    while (!to_exit)
    {
        /* Poll for and process events */
        glfwPollEvents();
        /* Render here */
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, out_image);
        glBlitFramebuffer(0, 0, nx, ny, 0, 0, nx, ny,
            GL_COLOR_BUFFER_BIT, GL_NEAREST);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        if (future.wait_for(std::chrono::seconds(0)) == std::future_status::ready || glfwWindowShouldClose(window))
        {
            to_exit = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
}

int main(int argc, const char **argv)
{
    constexpr int nx = 1024;
    constexpr int ny = 768;
    int ns = 100;
    constexpr int comp = 3; //RGB
    auto out_image = std::make_unique<GLubyte[]>(nx * ny * comp + 64);
    auto fout_image = std::make_unique<GLfloat[]>(nx * ny * comp + 64);
    memset(out_image.get(), 0, nx * ny * comp + 64);
    memset(fout_image.get(), 0.0f, nx * ny * comp + 64);
    //out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);
    bool to_exit = false;

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

    // Use cpp-taskflow https://github.com/cpp-taskflow/cpp-taskflow
    tf::Taskflow tf;

    // Start performance timer
    std::chrono::high_resolution_clock::time_point t11 = std::chrono::high_resolution_clock::now();

    // Initialize scene
    //std::unique_ptr<hitable> world(cornell_box_obj(cam, float(nx) / float(ny), lights));
    std::unique_ptr<Scene> scene(furnace_test_scene(float(nx) / float(ny)));

    std::chrono::high_resolution_clock::time_point t22 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_spann = std::chrono::duration_cast<std::chrono::duration<double>>(t22 - t11);
    std::cout << "\nBVH construction took me " << time_spann.count() << " seconds.";

    GLFWwindow* window;

    /* Register GLFW error callback */
    glfwSetErrorCallback(glfw_error_callback);
    /* Initialize GLFW */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    /* Also specify the OpenGL version (seems like this is sensitive to many potential issues) */
#if defined(__APPLE__) || defined(_WIN32)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
#else
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GL_FALSE);
#endif

#ifdef __APPLE__
    /* OSX requires forward compatibility for some reason */
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    window = glfwCreateWindow(nx, ny, "jam's ray tracer", NULL, NULL);
    if (!window) {
        fprintf(stderr, "ERROR: cannot open window with GLFW3\n");
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

#if defined(_WIN32) || defined(__linux__)
    /* Initialize GLEW */
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        // Problem: glewInit failed, something is seriously wrong.
        std::cout << "glewInit failed: " << glewGetErrorString(err) << std::endl;
        exit(1);
    }
#endif

    /* Create texture used to represent the color buffer */
    GLuint render_texture;
    glGenTextures(1, &render_texture);
    glBindTexture(GL_TEXTURE_2D, render_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);

    /* Create FBO to represent the framebuffer for blitting to the screen */
    GLuint fbo;
    glGenFramebuffers(1, &fbo);

    /* Bind FBO to both GL_FRAMEBUFFER and GL_READ_FRAMEBUFFER */
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, render_texture, 0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    /* Clear the window */
    glClear(GL_COLOR_BUFFER_BIT);

    // Use the renderer specified in template parameter
    path_prt renderer(scene.get(), ns);

    tf.parallel_for(ny - 1, 0, -1, [&](int j)
    {
        for (int i = 0; i < nx; i++)
        {
            //if (i <= 512)
            {
                if (to_exit) break;
                Vector3f col(0.0f, 0.0f, 0.0f);
                for (int s = 0; s < ns; s++)
                {
                    float u = float(i + gen_cano_rand()) / float(nx);
                    float v = float(j + gen_cano_rand()) / float(ny);
                    ray r = scene->cam.get_ray(u, v);
                    hit_record hrec;
                    const Vector3f sample = renderer.Li(r, scene.get(), 0, hrec, 0.0f);
                    assert(std::isfinite(sample[0])
                           && std::isfinite(sample[1])
                           && std::isfinite(sample[2]));
                    col += de_nan(sample);
                }
                col /= float(ns);

                float fr = col[0];
                float fg = col[1];
                float fb = col[2];

                int ir = std::min(int(pow(fr, 1.0f/2.2f) * 255.9999), 255);
                int ig = std::min(int(pow(fg, 1.0f/2.2f) * 255.9999), 255);
                int ib = std::min(int(pow(fb, 1.0f/2.2f) * 255.9999), 255);
                int index = (j * nx + i) * comp;

                // Store output pixels
                out_image[index] = (GLubyte)ir;
                out_image[index + 1] = (GLubyte)ig;
                out_image[index + 2] = (GLubyte)ib;

                fout_image[index] = fr;
                fout_image[index + 1] = fg;
                fout_image[index + 2] = fb;
            }
        }
        std::cout << ".";
    });

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    auto future = tf.dispatch();

    // Refresh window in background
    background_thread(future, out_image.get(), window, nx, ny, to_exit);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    std::cout << "\nIt took me " << time_span.count() << " seconds to render." << std::endl;

    std::cout << "Saving BMP..." << std::endl;
    image(out_image.get(), nx, ny, comp).save_image(formats::STBI_BMP);
    std::cout << "Saving JPG..." << std::endl;
    image(out_image.get(), nx, ny, comp).save_image(formats::STBI_JPG);
    std::cout << "Saving PNG..." << std::endl;
    image(out_image.get(), nx, ny, comp).save_image(formats::STBI_PNG);
    std::cout << "Saving PFM..." << std::endl;
    image_pfm(fout_image.get(), nx, ny, comp).save_image("out_test.pfm");

    glfwTerminate();

    return 0;
}
