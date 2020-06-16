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
    double Render(Scene *scene, viewer& film_viewer)
    {
        // Initialize viewer and create gl window
        GLFWwindow* window = film_viewer.init();
        // Use cpp-taskflow https://github.com/cpp-taskflow/cpp-taskflow
        tf::Taskflow tf;

        integrator::Render(scene, &film_viewer, tf);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        auto future = tf.dispatch();

        // Refresh window in background
        if constexpr (integrator::using_custom_viewer)
            integrator::background_thread(future, window, &film_viewer);
        else
            film_viewer.background_thread(future, window);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

        std::cout << "\nIt took me " << time_span.count() << " seconds to render." << std::endl;
        std::cout << "\n Number of ray traversals: " << Scene::num_rays_traversed;
        std::cout << "\n Number of ray-box intersections: " << Scene::num_ray_box_intersections;
        std::cout << "\n Number of ray-triangle intersections: " << Scene::num_ray_tri_intersections;
        std::cout << "\n Number of BVH nodes: " << Scene::num_bvh_nodes;
        std::cout << "\n Number of BVH leaves: " << Scene::num_bvh_leaf_nodes;
        std::cout << "\n Rays/second: " << (Scene::num_rays_traversed / time_span.count())/1e+6 << "M/sec" << std::endl;

        return time_span.count();
    }
};

#endif
