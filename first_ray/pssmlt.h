#ifndef PSSMLT_H
#define PSSMLT_H

#include "integrator.h"

/* Implements the PSSMLT paper "Simple and Robust Mutation Strategy for Metropolis Light Transport Algorithm"
   - Csaba Kelemen and László Szirmay-Kalos */
/* Reference: smallpssmlt by Professor Toshiya Hachisuka https://www.ci.i.u-tokyo.ac.jp/~hachisuka/smallpssmlt.cpp */
constexpr static int PixelWidth = 320;
constexpr static int PixelHeight = 180;
constexpr static int MaxPathLength = 20;
constexpr static int N_Init = 10000;
constexpr static float LargeStepProb = 0.3f;

constexpr static int NumRNGsPerEvent = 3;
constexpr static int MaxEvents = (MaxPathLength + 1);
constexpr static int NumStatesSubpath = ((MaxEvents + 2) * NumRNGsPerEvent);
constexpr static int NumStates = (NumStatesSubpath * 5);


// path data
struct Vert 
{ 
    Vector3f p, n; 
    hitable *hit_obj;
    Vert() {}; 
    Vert(Vector3f p_, Vector3f n_, hitable *hit_obj_) : p(p_), n(n_), hit_obj(hit_obj_) {}
};
struct PathContribution
{
    //x,y coordinate of eye ray
    float x, y;
    Vector3f c;
    //scalar contribution of path
    double sc;
    PathContribution() {};
    PathContribution(float x_, float y_, const Vector3f &c_, const double sc_) : x(x_), y(y_), c(c_), sc(sc_)
    {} 
};
struct Path
{
    Vert x[MaxEvents];
    PathContribution contrib;
    Vector3f camera_ray;
    int n;
    Path() { n = 0; }
};

struct pssmlt
{
    void TracePath(Path &path, const ray &r, Scene* scene);

    Path GenerateEyePath(const int MaxEyeEvents, Scene* scene);

    Vector3f Li(Path &path, const ray &r, Scene *scene, const int depth, const hit_record &prev_hrec,
        const float prev_bsdf_pdf);

    void Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf);

    constexpr static bool using_custom_viewer = true;

    void background_thread(const std::shared_future<void> &future, GLFWwindow *window, viewer *film_viewer);
};

#endif