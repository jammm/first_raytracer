#ifndef SPLIT_BVH_H
#define SPLIT_BVH_H

#include "hitable.h"

inline int box_x_compare(const void* a, const void* b)
{
	aabb box_left;
	aabb box_right;

	hitable* ah = *(hitable * *)a;
	hitable* bh = *(hitable * *)b;

	if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
		std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

	if (box_left.min.x() - box_right.min.x() < 0.0f)
		return -1;
	else
		return 1;
}

inline int box_y_compare(const void* a, const void* b)
{
	aabb box_left;
	aabb box_right;

	hitable* ah = *(hitable * *)a;
	hitable* bh = *(hitable * *)b;

	if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
		std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

	if (box_left.min.y() - box_right.min.y() < 0.0f)
		return -1;
	else
		return 1;
}

inline int box_z_compare(const void* a, const void* b)
{
	aabb box_left;
	aabb box_right;

	hitable* ah = *(hitable * *)a;
	hitable* bh = *(hitable * *)b;

	if (!ah->bounding_box(0.0f, 0.0f, box_left) || !bh->bounding_box(0.0f, 0.0f, box_right))
		std::cerr << "no bounding box in parallel_bvh_node constructor!\n";

	if (box_left.min.z() - box_right.min.z() < 0.0f)
		return -1;
	else
		return 1;
}

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

bool bvh_node::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    // Martin Lamber's BVH improvement https://twitter.com/Peter_shirley/status/1105292977423900673
    if (box.hit(r, t_min, t_max))
    {
        if (left->hit(r, t_min, t_max, rec))
        {
            right->hit(r, t_min, rec.t, rec);
            return true;
        }
        else
        {
            return right->hit(r, t_min, t_max, rec);
        }
    }
    return false;
}

// Construct BVH using SAH method
bvh_node::bvh_node(hitable **l, int n, float time0, float time1)
{
	std::unique_ptr<aabb[]> boxes(new aabb[n]);
	std::unique_ptr<float[]> left_area(new float[n]);
	std::unique_ptr<float[]> right_area(new float[n]);
	aabb main_box;
	l[0]->bounding_box(time0, time1, main_box);
	
	for (int i = 1; i < n; ++i)
	{
		aabb new_box;
		l[i]->bounding_box(time0, time1, new_box);
		main_box = surrounding_box(new_box, main_box);
	}

	int axis = main_box.longest_axis();
    if (axis == 0)
        qsort(l, n, sizeof(hitable *), box_x_compare);
    else if(axis == 1)
        qsort(l, n, sizeof(hitable *), box_y_compare);
    else
        qsort(l, n, sizeof(hitable *), box_z_compare);

	for (int i = 0; i < n; ++i)
		l[i]->bounding_box(time0, time1, boxes[i]);

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
		float SAH = i*left_area[i] + (n-i -1)*right_area[i+1];
		if (SAH < min_SAH)
		{
			min_SAH_idx = i;
			min_SAH = SAH;
		}
	}

    if (n == 1) {
        left = right = l[0];
    }
    else if (n == 2) {
        left = l[0];
        right = l[1];
    }
    else {
        left = new bvh_node(l, n / 2, time0, time1);
        right = new bvh_node(l + n / 2, n - n / 2, time0, time1);
    }

	box = main_box;
}

#endif