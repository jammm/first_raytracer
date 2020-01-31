#include "hitable_list.h"
#include "util.h"

bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const
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

bool hitable_list::bounding_box(float t0, float t1, aabb& box) const
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
		if (list[i]->bounding_box(t0, t1, temp_box))
		{
			box = surrounding_box(box, temp_box);
		}
		else
			return false;
	}

	return true;
}

float hitable_list::pdf_direct_sampling(const hit_record& lrec, const Vector3f& to_light) const
{
	float weight = 1.0 / list_size;
	float sum = 0;
	for (int i = 0; i < list_size; i++)
		sum += weight * list[i]->pdf_direct_sampling(lrec, to_light);

	return sum;
}

Vector3f hitable_list::sample_direct(hit_record& rec, const Vector3f& o, const Vector2f &sample) const
{
	float rand = sample[0];
	int index = int(rand * list_size);
	if (index == list_size) index -= 1;
	return list[index]->sample_direct(rec, o, sample);
}

int hitable_list::pick_sample(const float& sample) const
{
	int index = int(sample * list_size);
	if (index == list_size) index -= 1;

	return index;
}