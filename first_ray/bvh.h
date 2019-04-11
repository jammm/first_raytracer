#ifndef BVH_H
#define BVH_H

#include "hitable.h"

class bvh_node : public hitable
{
public:
	bvh_node() {}
	bvh_node(hitable **l, int n, float time0, float time1);
	virtual bool hit(const ray &r, float tmin, float tmax, hit_record &rec) const;
	virtual bool bounding_box(float t0, float t1, aabb &b) const;
	
	hitable *left;
	hitable *right;
	aabb box;
};

bool bvh_node::bounding_box(float t0, float t1, aabb &b) const
{
	b = box;
	return true;
}

#endif