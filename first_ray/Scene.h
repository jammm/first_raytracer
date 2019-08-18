#ifndef SCENE_H
#define SCENE_H

#include "camera.h"
#include "hitable_list.h"
#include "material.h"
#include <vector>

struct Scene
{
	std::unique_ptr<hitable> world;
	std::unique_ptr<environment_map> env_map;
	hitable_list lights;
	camera cam;

	Scene(hitable* world, environment_map* env_map, const camera& cam, const std::vector<hitable*> &lights) : world(world), env_map(env_map), cam(cam), lights(lights, lights.size())
	{}
};

#endif