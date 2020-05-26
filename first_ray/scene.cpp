#include "Scene.h"

std::atomic<unsigned int> Scene::num_rays_traversed = 0;
std::atomic<unsigned int> Scene::num_ray_tri_intersections = 0;
std::atomic<unsigned int> Scene::num_ray_box_intersections = 0;
std::atomic<unsigned int> Scene::num_bvh_nodes = 0;
std::atomic<unsigned int> Scene::num_bvh_leaf_nodes = 0;