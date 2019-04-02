
#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include "triangle.h"
#include "material.h"
#include "camera.h"
#include "image.h"
#include "texture.h"
#include "util.h"
#include <float.h>
#include <omp.h>
#include <thread>
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

void glfw_error_callback(int, const char* err_str)
{
    std::cout << "GLFW Error: " << err_str << std::endl;
}


hitable *random_scene()
{
    /* n == number of spheres */
    int n = 500;
    hitable **list = new hitable*[n + 1];
    //Large sphere's texture can be checkered
    texture *checker = new checker_texture(new constant_texture(Vector3f(0.2f, 0.3f, 0.1f)), new constant_texture(Vector3f(0.99f, 0.99f, 0.99f)));
    // TODO: Avoid hardcoded texture filenames. Use assimp!
    std::unique_ptr<image> img = std::make_unique<image>("cube/default.png");
    texture *image_tex = new image_texture(img);

    list[0] = new sphere(Vector3f(0, -1001, 0), 1000, new lambertian(checker));
    int i = 1;
    for (int a = -11;a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = drand48();
            Vector3f center(a + 0.9f * drand48(), 0.2f, b + 0.9f * drand48());
            if ((center - Vector3f(4, 0.2f, 0)).length() > 0.9f)
            {
                if (choose_mat < 0.8)
                {
                    //diffuse
                    list[i++] = new sphere(center, 0.2f, new lambertian(new constant_texture(Vector3f(drand48() * drand48(), drand48() * drand48(), drand48() * drand48()))));
                }
                else if (choose_mat < 0.95)
                {
                    //metal
                    list[i++] = new sphere(center, 0.2f,
                        new metal(Vector3f(0.5f*(1 + drand48()), 0.5f*(1 + drand48()), 0.5f*(1 + drand48())), 0.5f*drand48()));
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
    static std::vector <std::shared_ptr<hitable>> mesh = create_triangle_mesh("cube/cube.obj",
        std::make_shared<lambertian>(lambertian(image_tex)));

    for (auto triangle : mesh)
    {
        list[i++] = triangle.get();
    }

    //list[i++] = new sphere(Vector3f(0, 1, 0), 1.0f, new dielectric(1.5));
    //list[i++] = new sphere(Vector3f(4, 1, 0), 1.0f, new diffuse_light(new constant_texture(Vector3f(0.99f, 0.99f, 0.99f))));
    //list[i++] = new sphere(Vector3f(-4, 1, 0), 1.0f, new metal(Vector3f(0.7f, 0.6f, 0.5f), 0));
    //list[i++] = new triangle(Vector3f(0, 1, 0), Vector3f(4, 2, 0), Vector3f(-4, 1, 0), new diffuse_light(new constant_texture(Vector3f(0.99f, 0.99f, 0.99f))));

    return new hitable_list(list, i);
}

Vector3f color(const ray &r, hitable *world, int depth)
{
    hit_record rec;
    if (world->hit(r, 0.001f, FLT_MAX, rec))
    {
        ray scattered;
        Vector3f attenuation;
        Vector3f emitted = rec.mat_ptr->emitted(rec);
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        {
            return emitted + attenuation * color(scattered, world, depth + 1);
        }
        return emitted;
    }
    Vector3f unit_direction = unit_vector(r.direction());
    float t = 0.5f * (unit_direction.y() + 1.0f);
    return (1.0f - t)*Vector3f(1.0f, 1.0f, 1.0f) + t*Vector3f(0.5f, 0.7f, 1.0f);
    //return Vector3f(0, 0, 0);
}

int main()
{
    const int nx = 1024;
    const int ny = 768;
    const int ns = 100;
    const int comp = 3; //RGB
    GLubyte *out_image = new unsigned char[nx * ny * comp + 64];
    memset(out_image, 0, nx * ny * comp + 64);
    out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);
    bool to_exit = false;

    Vector3f lookfrom(13, 2, 3);
    Vector3f lookat(0, 0, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.05f;
    camera cam(lookfrom, lookat, Vector3f(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus);
    float R = (float) cos(M_PI / 4);
    int finished_threads = 0;

    // TODO: Read obj files for meshes. Scenes come later
    //hitable *list[5];
    //list[0] = new sphere(Vector3f(0.0f, 0.0f, -1.0f), 0.5f, new lambertian(Vector3f(0.1f, 0.2f, 0.5f)));
    //list[1] = new sphere(Vector3f(0.0f, -100.5f, -1.0f), 100.0f, new lambertian(Vector3f(0.8f, 0.8f, 0.0f)));
    //list[2] = new sphere(Vector3f(1.0f, 0, -1.0f), 0.5f, new metal(Vector3f(0.8f, 0.6f, 0.2f), 0.5f));
    //list[3] = new sphere(Vector3f(-1.0f, 0.0f, -1.0f), 0.5f, new dielectric(1.5f));
    //list[4] = new sphere(Vector3f(-1.0f, 0.0f, -1.0f), -0.45f, new dielectric(1.5f));
    //hitable *world = new hitable_list(list, 5);
    hitable *world = random_scene();

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
    std::cout << glGetError() << std::endl;
    glBindTexture(GL_TEXTURE_2D, render_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
    std::cout << glGetError() << std::endl;

    /* Create FBO to represent the framebuffer for blitting to the screen */
    GLuint fbo;
    glGenFramebuffers(1, &fbo);

    /* Bind FBO to both GL_FRAMEBUFFER and GL_READ_FRAMEBUFFER */
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, render_texture, 0);
    std::cout << glGetError() << std::endl;
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    // Flip output image vertically
    stbi_flip_vertically_on_write(true);

    /* Clear the window */
    glClear(GL_COLOR_BUFFER_BIT);

    // TODO: Parallelize this stuff
    // Use C++ std::thread, TBB or https://github.com/dougbinks/enkiTS

    #pragma omp parallel shared(to_exit)
    {
        #pragma omp for nowait
        for (int j = ny-1; j >= 0; j--)
        {
            for (int i = 0; i < nx; i++)
            {
                if (to_exit) break;
                Vector3f col(0.0f, 0.0f, 0.0f);
                for (int s = 0; s < ns; s++)
                {
                    float u = float(i + float(rand()) / float(RAND_MAX)) / float(nx);
                    float v = float(j + float(rand()) / float(RAND_MAX)) / float(ny);
                    ray r = cam.get_ray(u, v);
                    //Vector3f p = r.point_at_parameter(2.0f);
                    col += color(r, world, 0);
                }
                col /= float(ns);
                int ir = int(sqrt(col[0]) * 255.99);
                int ig = int(sqrt(col[1]) * 255.99);
                int ib = int(sqrt(col[2]) * 255.99);
                int index = (j * nx + i) * comp;

                // Store output pixels
                out_image[index]     = (unsigned char)ir;
                out_image[index + 1] = (unsigned char)ig;
                out_image[index + 2] = (unsigned char)ib;

            }
            if (omp_get_thread_num() == 0)
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
            }
            std::cout << ".";
            //std::cout << (float(ny - 1 - j) / float(ny)) * 100.0f << "%\r\r\r\r";
        }
        std::cout << finished_threads << std::endl;
        if(++finished_threads == omp_get_max_threads())
            to_exit = true;
        #pragma omp master
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
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
    }

    stbi_write_bmp("out_test.bmp", nx, ny, comp, (void *)out_image);
    stbi_write_png("out_test.png", nx, ny, comp, (void *)out_image, nx * comp);
    stbi_write_jpg("out_test.jpg", nx, ny, comp, (void *)out_image, 100);

    glfwTerminate();

    //TODO: Use unique_ptr to promote laziness
    delete[] out_image;
    delete[] world;

    return 0;
}
