#ifndef HITABLE_LIST_H_
#define HITABLE_LIST_H_

#include "hitable.h"

aabb surrounding_box(aabb box0, aabb box1)
{
	Vector3f small(fmin(box0.min.x(), box1.min.x()),
				   fmin(box0.min.y(), box1.min.y()),
				   fmin(box0.min.z(), box1.min.z()));
	Vector3f big  (fmin(box0.max.x(), box1.max.x()),
				   fmin(box0.max.y(), box1.max.y()),
				   fmin(box0.max.z(), box1.max.z()));

	return aabb(small, big);
}

class hitable_list : public hitable
{
public:
    hitable_list() {}
    hitable_list(hitable **l, int n) { list = l; list_size = n; }

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;
    hitable **list;
    int list_size;
};

bool hitable_list::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (int i = 0; i < list_size; i++)
    {
        if (list[i]->hit(r, t_min, (float)closest_so_far, temp_rec))
        {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
};

bool hitable_list::bounding_box(float t0, float t1, aabb &box) const
{
	if (list_size < 1) return false;
	aabb temp_box;
	bool first_true = list[0]->bounding_box(t0, t1, temp_box);
	if (first_true)
		return false;
	else
		box = temp_box;

	for (int i = 1; i < list_size; i++)
	{
		if (list[0]->bounding_box(t0, t1, temp_box))
		{
			box = surrounding_box(box, temp_box);
		}
		else
			return false;
	}

	return true;
}

#endif