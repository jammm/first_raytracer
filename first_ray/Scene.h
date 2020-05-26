#ifndef SCENE_H
#define SCENE_H

#include "camera.h"
#include "hitable_list.h"
#include "material.h"
#include <atomic>
#include <vector>

struct Scene
{
	std::unique_ptr<hitable> world;
	std::unique_ptr<environment_map> env_map;
	hitable_list lights;
	camera cam;

	// Add some metrics for profiling performance
	static std::atomic<unsigned int> num_rays_traversed;
	static std::atomic<unsigned int> num_ray_tri_intersections;
	static std::atomic<unsigned int> num_ray_box_intersections;
	static std::atomic<unsigned int> num_bvh_nodes;
	static std::atomic<unsigned int> num_bvh_leaf_nodes;

	Scene(hitable* world, environment_map* env_map, const camera& cam, const std::vector<hitable*> &lights) : world(world), env_map(env_map), lights(lights, lights.size()), cam(cam)
	{}
};

#endif