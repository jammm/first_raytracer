#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <cmath>
#include <array>
#include <memory>
#include "geometry.h"

constexpr float EPSILON = 1e-4f;
constexpr float SHADOW_EPSILON = 1e-3f;

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

inline float safe_sqrt(const float value)
{
    return std::max(0.0f, std::sqrt(value));
}

enum class microfacet_distributions
{
    ggx = 0,
    beckmann = 1
};

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

/* All these functions below are taken from mitsuba */
inline float modulo(float a, float b) {
    float r = std::fmod(a, b);
    return (r < 0.0f) ? r+b : r;
}

#if defined(_GNU_SOURCE)
inline void sincos(float theta, float *sin, float *cos) {
    sincosf(theta, sin, cos);
}

inline void sincos(double theta, double *sin, double *cos) {
    sincos(theta, sin, cos);
}

#else
inline void sincos(float theta, float* _sin, float *_cos) {
    *_sin = sinf(theta);
    *_cos = cosf(theta);
}

inline void sincos(double theta, double *_sin, double *_cos) {
    *_sin = sin(theta);
    *_cos = cos(theta);
}
#endif

#if defined(__LINUX__) && defined(__x86_64__)
inline float fastexp(float value) {
    return (float) exp((double) value);
}

inline double fastexp(double value) {
    return exp(value);
}

inline float fastlog(float value) {
    return (float) log((double) value);
}

inline double fastlog(double value) {
    return log(value);
}
#else
inline float fastexp(float value) {
    return expf(value);
}

inline double fastexp(double value) {
    return exp(value);
}

inline float fastlog(float value) {
    return logf(value);
}

inline double fastlog(double value) {
    return log(value);
}
#endif

inline float erfinv(float x) 
{
    // Based on "Approximating the erfinv function" by Mark Giles
    float w = -fastlog(((float)1 - x) * ((float)1 + x));
    float p;
    if (w < (float)5) 
    {
        w = w - (float)2.5;
        p = (float)2.81022636e-08;
        p = (float)3.43273939e-07 + p * w;
        p = (float)-3.5233877e-06 + p * w;
        p = (float)-4.39150654e-06 + p * w;
        p = (float)0.00021858087 + p * w;
        p = (float)-0.00125372503 + p * w;
        p = (float)-0.00417768164 + p * w;
        p = (float)0.246640727 + p * w;
        p = (float)1.50140941 + p * w;
    }
    else
    {
        w = std::sqrt(w) - (float)3;
        p = (float)-0.000200214257;
        p = (float)0.000100950558 + p * w;
        p = (float)0.00134934322 + p * w;
        p = (float)-0.00367342844 + p * w;
        p = (float)0.00573950773 + p * w;
        p = (float)-0.0076224613 + p * w;
        p = (float)0.00943887047 + p * w;
        p = (float)1.00167406 + p * w;
        p = (float)2.83297682 + p * w;
    }

    return p * x;
}

inline float erf_(double x) {
    double a1 = (float)0.254829592;
    double a2 = (float)-0.284496736;
    double a3 = (float)1.421413741;
    double a4 = (float)-1.453152027;
    double a5 = (float)1.061405429;
    double p = (float)0.3275911;

    // Save the sign of x
    double sign = copysignf(1.0f, x);
    x = std::abs(x);

    // A&S formula 7.1.26
    double t = (float)1.0 / ((float)1.0 + p * x);
    double y = (float)1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * fastexp(-x * x);

    return sign * y;
}

inline float hypot2(float a, float b)
{
    float r;
    if (std::abs(a) > std::abs(b)) {
        r = b / a;
        r = std::abs(a) * std::sqrt(1.0f + r * r);
    }
    else if (b != 0.0f) {
        r = a / b;
        r = std::abs(b) * std::sqrt(1.0f + r * r);
    }
    else {
        r = 0.0f;
    }
    return r;
}

#endif

