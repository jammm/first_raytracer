#ifndef PARALLEL_BVH_H
#define PARALLEL_BVH_H

#include "hitable.h"
#include "util.h"
#include <taskflow/taskflow.hpp>

class parallel_bvh_node : public hitable
{
public:
	parallel_bvh_node() {}
    parallel_bvh_node(tf::SubflowBuilder &subflow, hitable **l, int n, float time0, float time1);
	virtual bool hit(const ray &r, float tmin, float tmax, hit_record &rec) const;
	virtual bool bounding_box(float t0, float t1, aabb &b) const;

    static parallel_bvh_node *create_bvh(hitable **l, int n, float time0, float time1);

	hitable *left;
	hitable *right;
	aabb box;
};

bool parallel_bvh_node::bounding_box(float t0, float t1, aabb &b) const
{
	b = box;
	return true;
}

bool parallel_bvh_node::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
	if (box.hit(r, t_min, t_max))
	{
		hit_record left_rec, right_rec;
		bool hit_left = left->hit(r, t_min, t_max, left_rec);
		bool hit_right = right->hit(r, t_min, t_max, right_rec);

		if (hit_left && hit_right)
		{
			if (left_rec.t < right_rec.t)
				rec = left_rec;
			else
				rec = right_rec;
			return true;
		}
		else if (hit_left)
		{
			rec = left_rec;
			return true;
		}
		else if (hit_right)
		{
			rec = right_rec;
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

// Construct parallel_bvh using SAH method
parallel_bvh_node::parallel_bvh_node(tf::SubflowBuilder &subflow, hitable **l, int n, float time0, float time1)
    : left(nullptr), right(nullptr)
{
	std::unique_ptr<aabb[]> boxes(new aabb[n]);
	std::unique_ptr<float[]> left_area(new float[n]);
	std::unique_ptr<float[]> right_area(new float[n]);
	aabb main_box;
	bool dummy = l[0]->bounding_box(time0, time1, main_box);

	for (unsigned int i = 1; i < n; ++i)
	{
		aabb new_box;
		bool dummy = l[i]->bounding_box(time0, time1, new_box);
		main_box = surrounding_box(new_box, main_box);
	}

	int axis = main_box.longest_axis();
	if (axis == 0)
		qsort(l, n, sizeof(hitable *), box_x_compare);
	else if (axis == 1)
		qsort(l, n, sizeof(hitable *), box_y_compare);
	else
		qsort(l, n, sizeof(hitable *), box_z_compare);

	for (int i = 0; i < n; ++i)
		bool dummy = l[i]->bounding_box(time0, time1, boxes[i]);

	left_area[0] = boxes[0].area();
	aabb left_box = boxes[0];
	for (int i = 1; i < n - 1; ++i)
	{
		left_box = surrounding_box(left_box, boxes[i]);
		left_area[i] = left_box.area();
	}
	right_area[n - 1] = boxes[n - 1].area();
	aabb right_box = boxes[n - 1];
	for (int i = n - 2; i > 0; --i)
	{
		right_box = surrounding_box(right_box, boxes[i]);
		right_area[i] = right_box.area();
	}
	float min_SAH = FLT_MAX;
	int min_SAH_idx;
	for (int i = 0; i < n - 1; ++i)
	{
		float SAH = i * left_area[i] + (n - i - 1)*right_area[i + 1];
		if (SAH < min_SAH)
		{
			min_SAH_idx = i;
			min_SAH = SAH;
		}
	}

	if (min_SAH_idx == 0)
		left = l[0];
    else
    {
        subflow.emplace([=](tf::SubflowBuilder& subflow2)
        {
            left = new parallel_bvh_node(subflow2, l, min_SAH_idx + 1, time0, time1);
        }).name("left");
    }
	if (min_SAH_idx == n - 2)
		right = l[min_SAH_idx + 1];
    else
    {
        subflow.emplace([=](tf::SubflowBuilder& subflow2) {
            right = new parallel_bvh_node(subflow2, l + min_SAH_idx + 1, n - min_SAH_idx - 1, time0, time1);
        }).name("right");
    }

    subflow.detach();

	box = main_box;
}

parallel_bvh_node *parallel_bvh_node::create_bvh(hitable **l, int n, float time0, float time1)
{
    tf::Taskflow tf;
    parallel_bvh_node *root;
    auto task = tf.emplace([&](tf::SubflowBuilder &subflow)
    {
        root = new parallel_bvh_node(subflow, l, n, time0, time1);
    }).name("root");

    tf.dispatch().get();
    //tf.dump_topologies(std::cout);
    
    return root;
}

#endif