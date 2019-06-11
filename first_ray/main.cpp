
#include <iostream>
#include "hitable_list.h"
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
#include <float.h>
#include <taskflow/taskflow.hpp>
#include <chrono>

#ifdef _WIN32
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

inline float miWeight(float pdf1, float pdf2)
{
    pdf1 *= pdf1;
    pdf2 *= pdf2;
    return pdf1 / (pdf1 + pdf2);
}

void glfw_error_callback(int, const char* err_str)
{
    std::cout << "GLFW Error: " << err_str << std::endl;
}

hitable *random_scene(camera &cam, const float &aspect, std::vector<hitable *> &lights)
{
    /* n == number of spheres */
    int n = 1100000;
    hitable **list = new hitable*[n + 1];
    //Large sphere's texture can be checkered
    texture *checker = new checker_texture(new constant_texture(Vector3f(0.2f, 0.3f, 0.1f)), new constant_texture(Vector3f(0.99f, 0.99f, 0.99f)));
    // TODO: Avoid hardcoded texture filenames. Use assimp!
    std::unique_ptr<image> img = std::make_unique<image>("cube/default.png");
    texture *image_tex = new image_texture(img);

    list[0] = new sphere(Vector3f(0, -1002, 0), 1000, new lambertian(checker));
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

    // Load mesh from obj file
    // TODO: Change list to std::shared_ptr
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("CornellBox/CornellBox-Empty-CO.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    Vector3f lookfrom(12, 2, 3);
    Vector3f lookat(0, 0, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.001f;
    float vfov = 40.0f;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
    //return new bvh_node(list, i, 0.0f, 0.0f);
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
    float dist_to_focus = 10.0f;
    float aperture = 0.0f;
    float vfov = 40.0f;

    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    //return new hitable_list(list, i);
    return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
}

hitable *cornell_box_obj(camera &cam, const float &aspect, std::vector<hitable *> &lights)
{
    hitable **list = new hitable*[300];
    int i = 0;
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("CornellBox/CornellBox-Original.obj", lights);

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    Vector3f lookfrom(0, 1, 3.9f);
    Vector3f lookat(0, 1, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.0f;
    float vfov = 40.0f;
    cam = camera(lookfrom, lookat, Vector3f(0, 1, 0), vfov, aspect, aperture, dist_to_focus);

    return new hitable_list(std::vector<hitable *>(list, list + i), i);
    //return new bvh_node(list, i, 0.0f, 0.0f);
    //return parallel_bvh_node::create_bvh(list, i, 0.0f, 0.0f);
}

// TODO
// Use solid angle measure to make reference images
// Then use the reference against other measures/methods
Vector3f color(const ray &r, hitable *world, hitable *light_shape, int depth)
{
    hit_record hrec;
    if (world->hit(r, 1e-5, FLT_MAX, hrec))
    {
        scatter_record srec;
        Vector3f Li(0, 0, 0);
        const float invPi = 1 / M_PI;
        // Start with checking if the object hit is an emitter
        Li += hrec.mat_ptr->emitted(r, hrec);
        //Vector3f emitted = Vector3f(0, 0, 0);
        
        if (depth <= 50 && hrec.mat_ptr->scatter(r, hrec, srec))
        {
            if (srec.is_specular)
            {
                return srec.attenuation*color(srec.specular_ray, world, light_shape, depth + 1);
            }
            else
            {
                hitable_pdf plight(light_shape, hrec.p);
                //Vector3f on_light = Vector3f(0.129167, 1.98, 0.01882);
                Vector3f to_light = plight.generate();
                //float distance_squared = to_light.squared_length();
                float dist_to_light = to_light.length();
                to_light.make_unit_vector();

                //float light_cosine = -dot(to_light, Vector3f(0, -1, 0));
                //light_cosine = abs(light_cosine);

                /* Direct light sampling */
                ray shadow_ray = ray(hrec.p, to_light);
                hit_record lrec;

                // Calculate surface BSDF * cos(theta)
                float light_cosine = abs(dot(shadow_ray.d, hrec.normal));
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(hrec) * light_cosine;
                if (world->hit(shadow_ray, 1e-5, dist_to_light - 1e-5, lrec))
                {
                    if (dynamic_cast<diffuse_light *>(lrec.mat_ptr) != nullptr)
                    {
                        const float surface_bsdf_pdf = srec.pdf_ptr->value(srec.pdf_ptr->generate());
                        const float light_pdf = light_shape->pdf_direct_sampling(shadow_ray.o, shadow_ray.d);
                        const float weight = miWeight(surface_bsdf_pdf, light_pdf);

                        Li += lrec.mat_ptr->emitted(shadow_ray, lrec) * surface_bsdf * weight;
                    }
                }

                /* Sample BSDF to generate next ray direction for indirect lighting */
                ray wo(hrec.p, srec.pdf_ptr->generate());

                // srec.attenuation = bsdf weight == throughput
                return Li + srec.attenuation * color(wo, world, light_shape, depth + 1);


                //float invPi = 1 / M_PI;
                //Vector3f BRDF = srec.attenuation * invPi;
                //ray scattered(hrec.p, srec.pdf_ptr->generate());

                //Vector3f a = emitted * light_cosine * surface_cosine * BRDF * Visibility / (distance_squared);
                //Vector3f a = li_intensity * light_cosine * surface_cosine * Visibility;



                //return a + srec.attenuation * color(scattered, world, light_shape, depth + 1);
                //return a;



                /*
                Vector3f on_light = Vector3f(-0.130704165, 1.97999990, 0.152881131);
                Vector3f to_light = on_light - hrec.p;
                float distance_squared = to_light.squared_length();
                to_light.make_unit_vector();
                float light_cosine = dot(to_light, hrec.normal);
                if (light_cosine < 0)
                    return emitted;

                float light_area = 2.04014993;
                float pdf_light_area = 1 / light_area;
                float pdf_jacobian = light_cosine / distance_squared;

                ray scattered = ray(hrec.p, to_light);
                return emitted + srec.attenuation * hrec.mat_ptr->eval_bsdf(r, hrec, scattered)
                    * color(scattered, world, light_shape, depth + 1) * pdf_jacobian / pdf_light_area;
                */
            }
        }
        return Li;
    }
    //Vector3f unit_direction = unit_vector(r.direction());
    //float t = 0.5f * (unit_direction.y() + 1.0f);
    //return (1.0f - t)*Vector3f(1.0f, 1.0f, 1.0f) + t*Vector3f(0.5f, 0.7f, 1.0f);
    return Vector3f(0, 0, 0);
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

        if (glfwWindowShouldClose(window))
        {
            to_exit = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    }
}

int main()
{
    const int nx = 1024;
    const int ny = 768;
    const int ns = 100;
    const int comp = 3; //RGB
    auto out_image = std::make_unique<GLubyte[]>(nx * ny * comp + 64);
    auto fout_image = std::make_unique<GLfloat[]>(nx * ny * comp + 64);
    memset(out_image.get(), 0, nx * ny * comp + 64);
    memset(fout_image.get(), 0.0f, nx * ny * comp + 64);
    //out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);
    bool to_exit = false;

    // Use cpp-taskflow https://github.com/cpp-taskflow/cpp-taskflow
    tf::Taskflow tf;

    // Initialize scene
    camera cam;
    std::vector<hitable *> lights;

    std::chrono::high_resolution_clock::time_point t11 = std::chrono::high_resolution_clock::now();
    std::unique_ptr<hitable> world(cornell_box_obj(cam, float(nx) / float(ny), lights));

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
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GL_FALSE);

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

#ifdef _WIN32
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

    // TODO: Find a better way to specify lights in the scene
    hitable_list hlist(lights, 2);

    tf.parallel_for(ny - 1, 0, -1, [&](int j)
    {
        for (int i = 0; i < nx; i++)
        {
            //if (i == 783 && j == 411)
            {
                if (to_exit) break;
                Vector3f col(0.0f, 0.0f, 0.0f);
                for (int s = 0; s < ns; s++)
                {
                    float u = float(i + gen_cano_rand()) / float(nx);
                    float v = float(j + gen_cano_rand()) / float(ny);
                    ray r = cam.get_ray(u, v);
                    col += de_nan(color(r, world.get(), &hlist, 0));
                }
                col /= float(ns);

                float fr = col[0];
                float fg = col[1];
                float fb = col[2];

                int ir = std::min(int(sqrt(fr) * 255.99), 255);
                int ig = std::min(int(sqrt(fg) * 255.99), 255);
                int ib = std::min(int(sqrt(fb) * 255.99), 255);
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
        //std::cout << (float(ny - 1 - j) / float(ny)) * 100.0f << "%\r\r\r\r";
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
