
#include <iostream>
#include "hitable_list.h"
#include "sphere.h"
#include "material.h"
#include "camera.h"
#include "texture.h"
#include "util.h"
#include <float.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"




hitable *random_scene()
{
    /* n == number of spheres */
    int n = 500;
    hitable **list = new hitable*[n + 1];
    //Large sphere's texture can be checkered
    texture *checker = new checker_texture(new constant_texture(vec3(0.2f, 0.3f, 0.1f)), new constant_texture(vec3(0.9f, 0.9f, 0.9f)));
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
    int i = 1;
    for (int a = -11;a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = drand48();
            vec3 center(a + 0.9f * drand48(), 0.2f, b + 0.9f * drand48());
            if ((center - vec3(4, 0.2f, 0)).length() > 0.9f)
            {
                if (choose_mat < 0.8)
                {
                    //diffuse
                    list[i++] = new sphere(center, 0.2f, new lambertian(new constant_texture(vec3(drand48() * drand48(), drand48() * drand48(), drand48() * drand48()))));
                }
                else if (choose_mat < 0.95)
                {
                    //metal
                    list[i++] = new sphere(center, 0.2f,
                        new metal(vec3(0.5f*(1 + drand48()), 0.5f*(1 + drand48()), 0.5f*(1 + drand48())), 0.5f*drand48()));
                }
                else
                {
                    //dielectric
                    list[i++] = new sphere(center, 0.2f, new dielectric(1.5f));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0f, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0f, new lambertian(new constant_texture(vec3(0.4f, 0.2f, 0.1f))));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0f, new metal(vec3(0.7f, 0.6f, 0.5f), 0));

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
    const int nx = 1024;
    const int ny = 768;
    const int ns = 100;
    const int comp = 3; //RGB
    GLubyte *out_image = new unsigned char[nx * ny * comp + 64];
    memset(out_image, 0, nx * ny * comp + 64);
    out_image = (GLubyte *)(((std::size_t)out_image) >> 6 <<6);
    bool to_exit = false;

    vec3 lookfrom(13, 2, 3);
    vec3 lookat(0, 0, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.1f;
    camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus);
    float R = (float) cos(M_PI / 4);

    // TODO: Read obj files for meshes. Scenes come later
    //hitable *list[5];
    //list[0] = new sphere(vec3(0.0f, 0.0f, -1.0f), 0.5f, new lambertian(vec3(0.1f, 0.2f, 0.5f)));
    //list[1] = new sphere(vec3(0.0f, -100.5f, -1.0f), 100.0f, new lambertian(vec3(0.8f, 0.8f, 0.0f)));
    //list[2] = new sphere(vec3(1.0f, 0, -1.0f), 0.5f, new metal(vec3(0.8f, 0.6f, 0.2f), 0.5f));
    //list[3] = new sphere(vec3(-1.0f, 0.0f, -1.0f), 0.5f, new dielectric(1.5f));
    //list[4] = new sphere(vec3(-1.0f, 0.0f, -1.0f), -0.45f, new dielectric(1.5f));
    //hitable *world = new hitable_list(list, 5);
    hitable *world = random_scene();

    GLFWwindow* window;

    /* Initialize GLFW */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow(nx, ny, "jam's ray tracer", NULL, NULL);
    if (!window) {
        fprintf(stderr, "ERROR: cannot not open window with GLFW3\n");
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    /* Initialize GLEW */
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        // Problem: glewInit failed, something is seriously wrong.
        std::cout << "glewInit failed: " << glewGetErrorString(err) << std::endl;
        exit(1);
    }

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
    // Use C++ std::thread or https://github.com/dougbinks/enkiTS
    for (int j = ny-1; j >= 0; j--)
    {
        if (to_exit) break;
        for (int i = 0; i < nx; i++)
        {
            if (to_exit) break;
            vec3 col(0.0f, 0.0f, 0.0f);
            for (int s = 0; s < ns; s++)
            {
                if (glfwWindowShouldClose(window))
                {
                    to_exit = true;
                    break;
                }

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

            // Store output pixels
            out_image[index]     = unsigned char(ir);
            out_image[index + 1] = unsigned char(ig);
            out_image[index + 2] = unsigned char(ib);

            /* Render here */
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, out_image);
            glBlitFramebuffer(0, 0, nx, ny, 0, 0, nx, ny,
                GL_COLOR_BUFFER_BIT, GL_NEAREST);

            /* Swap front and back buffers */
            glfwSwapBuffers(window);

            /* Poll for and process events */
            glfwPollEvents();
        }
        std::cout << ".";
        //std::cout << (float(ny - 1 - j) / float(ny)) * 100.0f << "%\r\r\r\r";
    }

    stbi_write_bmp("out_test.bmp", nx, ny, comp, (void *)out_image);
    stbi_write_png("out_test.png", nx, ny, comp, (void *)out_image, nx * comp);
    stbi_write_jpg("out_test.jpg", nx, ny, comp, (void *)out_image, 100);

    glfwTerminate();
    delete out_image;

    return 0;
}