#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <cmath>
#include <array>
#include <memory>
#include "geometry.h"

constexpr float EPSILON = 1e-5f;
constexpr float SHADOW_EPSILON = 1e-4f;

inline float unit_angle(const Vector3f& u, const Vector3f& v) {
	if (dot(u, v) < 0)
		return M_PI - 2 * std::asin(0.5f * (v + u).length());
	else
		return 2 * std::asin(0.5f * (v - u).length());
}

inline Vector2f random_in_unit_disk(Vector2f sample)
{
	Vector2f ab = sample * 2.0f - Vector2f(1.0f, 1.0f);

	if (ab[0] == 0.0f && ab[1] == 0.0f) { return Vector2f(0.0f, 0.0f); }

	Vector2f ab2 = ab * ab;
	float phi, r;
	if (ab2[0] > ab2[1]) { // use squares instead of absolute values
		r = ab[0];
		phi = (M_PI / 4.0f) * (ab[1] / ab[0]);
	}
	else {
		r = ab[1];
		phi = (M_PI / 2.0f) - (M_PI / 4.0f) * (ab[0] / ab[1]);
	}
	float cosphi = std::cos(phi);
	float sinphi = std::sin(phi);

	return r * Vector2f(cosphi, sinphi);
}

// MIS weight using power heuristic
inline float miWeight(float pdf1, float pdf2)
{
    pdf1 *= pdf1;
    pdf2 *= pdf2;
    return pdf1 / (pdf1 + pdf2);
}

static inline float FromSrgb(const float &v)
{
    if (v <= 0.04045) return v * (1.0 / 12.92);
    return std::pow((v + 0.055) * (1.0 / 1.055), 2.4);
}

static inline Vector3f FromSrgb(const Vector3f &v)
{
	return Vector3f(FromSrgb(v.e[0]), FromSrgb(v.e[1]), FromSrgb(v.e[2]));
}

inline Vector3f reflect(const Vector3f &v, const Vector3f &n)
{
    return v - 2 * dot(v, n)*n;
}

#endif

