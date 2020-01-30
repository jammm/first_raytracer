#ifndef GL_INCLUDES_H
#define GL_INCLUDES_H

#if defined(_WIN32) || defined(__linux__)
#define GLEW_STATIC
#include <GL/glew.h>
#define _CRT_SECURE_NO_DEPRECATE
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#endif
#include <GLFW/glfw3.h>

#endif