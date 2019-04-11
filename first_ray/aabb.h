#ifndef AABB_H
#define AABB_H

#include "ray.h"

struct aabb
{
public:
	aabb() {}
	aabb(const Vector3f &min, const Vector3f &max) : min(min), max(max) {}


	// AABB intersection using slab method as described by peter shirley's ray tracing the next week book.
	bool hit(const ray &r, float tmin, float tmax)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			float invD = 1.0f / r.direction()[i];
			float t0 = (min[i] - r.origin()[i]) * invD;
			float t1 = (max[o] - r.origin()[i]) * invD;

			if (t1 < t0)
				std::swap(t1, t0);
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
				return false;
		}
	}

	Vector3f min;
	Vector3f max;

};

#endif