#include "viewer.h"
#include "image.h"
#include <iostream>

static void glfw_error_callback(int, const char* err_str)
{
    std::cout << "GLFW Error: " << err_str << std::endl;
}

viewer::viewer(const int nx, const int ny, const uint_fast64_t ns, const int num_channels) : to_exit(false), nx(nx), ny(ny), ns(ns), num_channels(num_channels)
{
    out_image = std::make_unique<GLubyte[]>(nx * ny * num_channels + 64);
    fout_image = std::make_unique<GLdouble[]>(nx * ny * num_channels + 64);
    memset(out_image.get(), 0, nx * ny * num_channels + 64);
    memset(fout_image.get(), 0.0, nx * ny * num_channels + 64);
}

GLFWwindow* viewer::init()
{
    /* Register GLFW error callback */
    glfwSetErrorCallback(glfw_error_callback);
    /* Initialize GLFW */
    if (!glfwInit())
        return nullptr;

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
        return nullptr;
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

    return window;
}

void viewer::background_thread(const std::shared_future<void>& future, GLFWwindow* window)
{
    while (!to_exit)
    {
        /* Poll for and process events */
        glfwPollEvents();
        /* Render here */
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB, GL_UNSIGNED_BYTE, out_image.get());
        glBlitFramebuffer(0, 0, nx, ny, 0, 0, nx, ny,
            GL_COLOR_BUFFER_BIT, GL_NEAREST);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        if (future.wait_for(std::chrono::seconds(0)) == std::future_status::ready || glfwWindowShouldClose(window))
        {
            to_exit = true;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(9));
    }
}

void viewer::add_sample(const Vector2i& pixel, Vector3f &sample)
{
    sample /= double(ns);

    const double &fr = sample[0];
    const double &fg = sample[1];
    const double &fb = sample[2];

    const int ir = std::min(int(pow(fr, 1.0 / 2.2f) * 255.9999), 255);
    const int ig = std::min(int(pow(fg, 1.0 / 2.2f) * 255.9999), 255);
    const int ib = std::min(int(pow(fb, 1.0 / 2.2f) * 255.9999), 255);
    const int index = (pixel[1] * nx + pixel[0]) * num_channels;

    // Store output pixels
    out_image[index] = (GLubyte)ir;
    out_image[index + 1] = (GLubyte)ig;
    out_image[index + 2] = (GLubyte)ib;

    fout_image[index] = fr;
    fout_image[index + 1] = fg;
    fout_image[index + 2] = fb;
}

void viewer::save_and_destroy(const std::string &filename)
{
    std::cout << "Saving BMP..." << std::endl;
    image(out_image.get(), nx, ny, num_channels).save_image(filename, formats::STBI_BMP);
    std::cout << "Saving JPG..." << std::endl;
    image(out_image.get(), nx, ny, num_channels).save_image(filename, formats::STBI_JPG);
    std::cout << "Saving PNG..." << std::endl;
    image(out_image.get(), nx, ny, num_channels).save_image(filename, formats::STBI_PNG);
    std::cout << "Saving PFM..." << std::endl;
    image_pfm(fout_image.get(), nx, ny, num_channels).save_image(filename + ".pfm");

    glfwTerminate();
}
