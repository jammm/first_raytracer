#ifndef UTIL_H
#define UTIL_H

#include <random>
#include <cmath>
#include <array>
#include <memory>
#include "geometry.h"

constexpr double EPSILON = 1e-4;
constexpr double SHADOW_EPSILON = 1e-3f;
constexpr double DELTA_EPSILON = 1e-3f;

inline double unit_angle(const Vector3f& u, const Vector3f& v) {
	if (dot(u, v) < 0)
		return M_PI - 2 * std::asin(0.5 * (v + u).length());
	else
		return 2 * std::asin(0.5 * (v - u).length());
}

inline Vector2f random_in_unit_disk(Vector2f sample)
{
	Vector2f ab = sample * 2.0 - Vector2f(1.0, 1.0);

	if (ab[0] == 0.0 && ab[1] == 0.0) { return Vector2f(0.0, 0.0); }

	Vector2f ab2 = ab * ab;
	double phi, r;
	if (ab2[0] > ab2[1]) { // use squares instead of absolute values
		r = ab[0];
		phi = (M_PI / 4.0) * (ab[1] / ab[0]);
	}
	else {
		r = ab[1];
		phi = (M_PI / 2.0) - (M_PI / 4.0) * (ab[0] / ab[1]);
	}
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);

	return r * Vector2f(cosphi, sinphi);
}

inline double safe_sqrt(const double value)
{
    return std::max(0.0, std::sqrt(value));
}

enum class microfacet_distributions
{
    ggx = 0,
    beckmann = 1
};

// MIS weight using power heuristic
inline double miWeight(double pdf1, double pdf2)
{
    pdf1 *= pdf1;
    pdf2 *= pdf2;
    return pdf1 / (pdf1 + pdf2);
}

static inline double FromSrgb(const double &v)
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
    return unit_vector(v - 2 * dot(v, n)*n);
}

/* All these functions below are taken from mitsuba */
inline Vector3f refract(const Vector3f &wi, const Vector3f &n, double eta, double cosThetaT) {
    if (cosThetaT < 0)
        eta = 1 / eta;

    return unit_vector(n * (dot(wi, n) * eta + cosThetaT) - wi * eta);
}

inline double fresnelDielectricExt(double cosThetaI_, double &cosThetaT_, double eta) {
    if (eta == 1) 
    {
        cosThetaT_ = -cosThetaI_;
        return 0.0;
    }

    /* Using Snell's law, calculate the squared sine of the
       angle between the normal and the transmitted ray */
    double scale = (cosThetaI_ > 0) ? 1/eta : eta,
          cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

    /* Check for total internal reflection */
    if (cosThetaTSqr <= 0.0)
    {
        cosThetaT_ = 0.0;
        return 1.0;
    }

    /* Find the absolute cosines of the incident/transmitted rays */
    double cosThetaI = std::abs(cosThetaI_);
    double cosThetaT = std::sqrt(cosThetaTSqr);

    double Rs = (cosThetaI - eta * cosThetaT)
             / (cosThetaI + eta * cosThetaT);
    double Rp = (eta * cosThetaI - cosThetaT)
             / (eta * cosThetaI + cosThetaT);

    cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

    /* No polarization -- return the unpolarized reflectance */
    return 0.5 * (Rs * Rs + Rp * Rp);
}

inline double modulo(double a, double b) {
    double r = std::fmod(a, b);
    return (r < 0.0) ? r+b : r;
}

inline int modulo(int a, int b) {
    int r = a % b;
    return (r < 0) ? r+b : r;
}

#if defined(_GNU_SOURCE)
inline void sincos(float theta, float *sin, float *cos) {
    sincosf(theta, sin, cos);
}

inline void sincos(double theta, double *sin, double *cos) {
    sincos(theta, sin, cos);
}

#else
inline void sincos(double theta, double* _sin, double *_cos) {
    *_sin = sinf(theta);
    *_cos = cosf(theta);
}

inline void sincos(float theta, float *_sin, float *_cos) {
    *_sin = sin(theta);
    *_cos = cos(theta);
}
#endif

#if defined(__LINUX__) && defined(__x86_64__)
inline float fastexp(float value) {
    return (float) exp((float) value);
}

inline double fastexp(double value) {
    return exp(value);
}

inline float fastlog(float value) {
    return (float) log((float) value);
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

inline double erfinv(double x) 
{
    // Based on "Approximating the erfinv function" by Mark Giles
    double w = -fastlog(((double)1 - x) * ((double)1 + x));
    double p;
    if (w < (double)5) 
    {
        w = w - (double)2.5;
        p = (double)2.81022636e-08;
        p = (double)3.43273939e-07 + p * w;
        p = (double)-3.5233877e-06 + p * w;
        p = (double)-4.39150654e-06 + p * w;
        p = (double)0.00021858087 + p * w;
        p = (double)-0.00125372503 + p * w;
        p = (double)-0.00417768164 + p * w;
        p = (double)0.246640727 + p * w;
        p = (double)1.50140941 + p * w;
    }
    else
    {
        w = std::sqrt(w) - (double)3;
        p = (double)-0.000200214257;
        p = (double)0.000100950558 + p * w;
        p = (double)0.00134934322 + p * w;
        p = (double)-0.00367342844 + p * w;
        p = (double)0.00573950773 + p * w;
        p = (double)-0.0076224613 + p * w;
        p = (double)0.00943887047 + p * w;
        p = (double)1.00167406 + p * w;
        p = (double)2.83297682 + p * w;
    }

    return p * x;
}

inline double erf_(double x) {
    double a1 = (double)0.254829592;
    double a2 = (double)-0.284496736;
    double a3 = (double)1.421413741;
    double a4 = (double)-1.453152027;
    double a5 = (double)1.061405429;
    double p = (double)0.3275911;

    // Save the sign of x
    double sign = copysignf(1.0, x);
    x = std::abs(x);

    // A&S formula 7.1.26
    double t = (double)1.0 / ((double)1.0 + p * x);
    double y = (double)1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * fastexp(-x * x);

    return sign * y;
}

inline double hypot2(double a, double b)
{
    double r;
    if (std::abs(a) > std::abs(b)) {
        r = b / a;
        r = std::abs(a) * std::sqrt(1.0 + r * r);
    }
    else if (b != 0.0) {
        r = a / b;
        r = std::abs(b) * std::sqrt(1.0 + r * r);
    }
    else {
        r = 0.0;
    }
    return r;
}

#endif

