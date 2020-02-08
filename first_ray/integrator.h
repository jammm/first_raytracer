#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "hitable_list.h"
#include "geometry.h"
#include "Scene.h"
#include "viewer.h"
#include "sampler.h"
#include <taskflow/taskflow.hpp>

template <typename integrator>
struct renderer : public integrator
{
    void Render(Scene *scene, viewer& film_viewer)
    {
        // Initialize viewer and create gl window
        GLFWwindow* window = film_viewer.init();
        // Use cpp-taskflow https://github.com/cpp-taskflow/cpp-taskflow
        tf::Taskflow tf(1);

        integrator::Render(scene, &film_viewer, tf);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        auto future = tf.dispatch();

        // Refresh window in background
        film_viewer.background_thread(future, window);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

        std::cout << "\nIt took me " << time_span.count() << " seconds to render." << std::endl;
    }
};

#endif
