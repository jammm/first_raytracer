#ifndef AABB_H
#define AABB_H

#include "ray.h"

class aabb
{
public:
	aabb() {}
	aabb(const Vector3f &min, const Vector3f &max) : min(min), max(max), size(max - min) {}


	// AABB intersection using slab method as described by peter shirley's ray tracing the next week book.
	bool hit(const ray &r, double tmin, double tmax) const
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			double invD = 1.0 / r.direction()[i];
			double t0 = (min[i] - r.origin()[i]) * invD;
			double t1 = (max[i] - r.origin()[i]) * invD;

			if (invD < 0.0)
				std::swap(t0, t1);
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax < tmin)
				return false;
		}
        
        return true;
	}

	int longest_axis() const
	{
		int axis = 0;

		if (size.y() > size.x())
			axis = 1;
		else if (size.z() > size.x())
			axis = 2;

		return axis;
	}

	// Returns surface area of AABB
	inline double area() const
	{
		return 2 * ((size.x() * size.y()) + (size.y() * size.z()) + (size.x() * size.z()));
	}

	Vector3f min;
	Vector3f max;
	Vector3f size;

};

inline aabb surrounding_box(aabb box0, aabb box1)
{
    Vector3f small(fmin(box0.min.x(), box1.min.x()),
                   fmin(box0.min.y(), box1.min.y()),
                   fmin(box0.min.z(), box1.min.z()));

    Vector3f big(  fmax(box0.max.x(), box1.max.x()),
                   fmax(box0.max.y(), box1.max.y()),
                   fmax(box0.max.z(), box1.max.z()));

    return aabb(small, big);
}

#endif
