#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <cmath>
#include <array>
#include <memory>
#include "geometry.h"

constexpr float EPSILON = 1e-4;
constexpr float SHADOW_EPSILON = 1e-3;

// Generates the canonical uniform random variable ξ
inline float gen_cano_rand()
{
    thread_local static std::random_device seed_gen;
    thread_local static std::mt19937 engine(seed_gen());
    thread_local static std::uniform_real_distribution<> dist(0.0f, 1.0f);
    return dist(engine);
}

inline float unit_angle(const Vector3f& u, const Vector3f& v) {
	if (dot(u, v) < 0)
		return M_PI - 2 * std::asin(0.5f * (v + u).length());
	else
		return 2 * std::asin(0.5f * (v - u).length());
}

inline Vector3f random_in_unit_disk()
{
	Vector3f p;
	do
	{
		p = 2.0f * Vector3f(gen_cano_rand(), gen_cano_rand(), 0.0f) - Vector3f(1.0f, 1.0f, 0.0f);
	} while (dot(p, p) >= 1.0f);

	return p;
}

// MIS weight using power heuristic
inline float miWeight(float pdf1, float pdf2)
{
    pdf1 *= pdf1;
    pdf2 *= pdf2;
    return pdf1 / (pdf1 + pdf2);
}

static inline float FromSrgb(const float v)
{
    if (v <= 0.04045) return v * (1.0 / 12.92);
    return std::pow((v + 0.055) * (1.0 / 1.055), 2.4);
}

inline Vector3f reflect(const Vector3f &v, const Vector3f &n)
{
    return v - 2 * dot(v, n)*n;
}

#endif

