#pragma once

#include "ray.h"

#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__ ((aligned (16)))
#endif

// turn those verbose intrinsics into something readable.
#define loadps(mem)			_mm_load_ps((const double * const)(mem))
#define storess(ss,mem)		_mm_store_ss((double * const)(mem),(ss))
#define minss				_mm_min_ss
#define maxss				_mm_max_ss
#define minps				_mm_min_ps
#define maxps				_mm_max_ps
#define mulps				_mm_mul_ps
#define subps				_mm_sub_ps
#define rotatelps(ps)		_mm_shuffle_ps((ps),(ps), 0x39)	// a,b,c,d -> b,c,d,a
#define muxhps(low,high)	_mm_movehl_ps((low),(high))	// low{a,b,c,d}|high{e,f,g,h} = {c,d,g,h}


static const double flt_plus_inf = -logf(0);	// let's keep C and C++ compilers happy.
static const double _MM_ALIGN16
	ps_cst_plus_inf[4] = { flt_plus_inf,  flt_plus_inf,  flt_plus_inf,  flt_plus_inf },
	ps_cst_minus_inf[4] = { -flt_plus_inf, -flt_plus_inf, -flt_plus_inf, -flt_plus_inf };

struct aabb
{
	aabb() {}
	aabb(const Vector3f &bmin, const Vector3f &bmax) : min(bmin), max(bmax), size(max - min) {}

	// AABB intersection using slab method as described by peter shirley's ray tracing the next week book.
	bool __vectorcall hit(const ray &r, double tmin, double tmax) const
	{
		// you may already have those values hanging around somewhere
		static const double _MM_ALIGN16 tmax_arr[4] = { tmax, tmax, tmax, tmax };
		static const double _MM_ALIGN16 tmin_arr[4] = { tmin, tmin, tmin, tmin };

		const __m128
			tmax4 = loadps(tmax_arr),
			tmin4 = loadps(tmin_arr);

		// use whatever's apropriate to load.
		Vector3f i_dir = Vector3f(1,1,1)/r.d;
		const __m128
			box_min = loadps(&min),
			box_max = loadps(&max),
			pos = loadps(&r.o),
			inv_dir = loadps(&i_dir);

		// use a div if inverted directions aren't available
		const __m128 l1 = mulps(subps(box_min, pos), inv_dir);
		const __m128 l2 = mulps(subps(box_max, pos), inv_dir);

		// the order we use for those min/max is vital to filter out
		// NaNs that happens when an inv_dir is +/- inf and
		// (box_min - pos) is 0. inf * 0 = NaN
		const __m128 filtered_l1a = minps(l1, tmax4);
		const __m128 filtered_l2a = minps(l2, tmax4);

		const __m128 filtered_l1b = maxps(l1, tmin4);
		const __m128 filtered_l2b = maxps(l2, tmin4);

		// now that we're back on our feet, test those slabs.
		__m128 lmax = maxps(filtered_l1a, filtered_l2a);
		__m128 lmin = minps(filtered_l1b, filtered_l2b);

		// unfold back. try to hide the latency of the shufps & co.
		const __m128 lmax0 = rotatelps(lmax);
		const __m128 lmin0 = rotatelps(lmin);
		lmax = minss(lmax, lmax0);
		lmin = maxss(lmin, lmin0);

		const __m128 lmax1 = muxhps(lmax, lmax);
		const __m128 lmin1 = muxhps(lmin, lmin);
		lmax = minss(lmax, lmax1);
		lmin = maxss(lmin, lmin1);

		const bool ret = _mm_comige_ss(lmax, _mm_setzero_ps()) & _mm_comige_ss(lmax, lmin);

		storess(lmin, &tmin);
		storess(lmax, &tmax);

		return ret;
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