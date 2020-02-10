#ifndef PDF_H_
#define PDF_H_

#include "geometry.h"
#include "onb.h"
#include "hitable.h"
#include "util.h"
#include "microfacet.h"
#include <algorithm>

class rough_conductor;

inline Vector3f hemisphere_to_cosine_direction(const Vector2f &sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float r = sqrt(r0);
    const float phi = 2 * (float)M_PI * r1;
    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}

inline Vector3f hemisphere_to_cosine_direction(float &theta, float& phi, const Vector2f &sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float r = sqrt(r0);
    theta = asin(r);
    phi = 2 * (float)M_PI * r1;
    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}


inline Vector3f hemisphere_to_cosine_power_direction(const Vector2f& sample)
{
    const float &r0 = sample[0];
    const float &r1 = sample[1];
    const float sin_theta = sqrt(1 - r0);
    const float r = 2 * (float)M_PI *r1;

    return Vector3f(sin_theta*cos(r), sin_theta * sin(r), sqrt(r1));
}

inline Vector3f uniform_sample_sphere(const Vector2f &sample) {
    Vector2f u(sample);
    const float  z = 1 - 2 * u[0];
    const float r = std::sqrt(std::max((float)0, (float)1 - z * z));
    const float phi = 2 * (float)M_PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline Vector3f random_to_sphere(const float &radius, const float &distance_squared, const Vector2f& sample)
{
    const float &r1 = sample[0];
    const float &r2 = sample[1];
    float z = 1 + r2 * (sqrt(1 - radius*radius/distance_squared) - 1);
    float phi = 2*(float)M_PI*r1;
    float x = cos(phi)*sqrt(1-z*z);
    float y = sin(phi)*sqrt(1-z*z);

    return Vector3f(x, y, z);
}

class pdf
{
public:
    virtual float value(const hit_record &, const Vector3f &) const = 0;
    virtual Vector3f generate(const Vector2f &sample, const hit_record &hrec) const = 0;
    virtual ~pdf() = 0;
};

class cosine_pdf : public pdf
{
public:
    cosine_pdf(const Vector3f &w) { uvw.build_from_w(w); }
    float value(const hit_record &hrec, const Vector3f &direction) const override
    {
        float cosine = std::max<float>(dot(uvw.w(), unit_vector(direction)), 0.0f);
        
        return cosine / (float) (float)M_PI;
    }

    Vector3f generate(const Vector2f &sample, const hit_record &hrec) const override
    {
        return uvw.fromLocal(hemisphere_to_cosine_direction(sample));
    }

    onb uvw;
};

class cosine_power_pdf : public pdf
{
public:
    cosine_power_pdf(const ray &r_in, const Vector3f &w, const float &specular_exponent) : specular_exponent(specular_exponent), r_in(r_in) { uvw.build_from_w(w); }
    // TODO: modify this stuff - modified phong BRDF paper
    float value(const hit_record& hrec, const Vector3f& wo) const override
    {
        if (dot(hrec.normal, wo) <= 0
            || dot(hrec.normal, hrec.wi) <= 0)
            return 0.0f;
        const double alpha = std::max(0.0f, dot(reflect(-hrec.wi, uvw.w()), wo));
        const double specular_pdf = std::pow(alpha, (double)specular_exponent);

        return specular_pdf * (specular_exponent + 1.0f) / (2 * (float)M_PI);
    }

    Vector3f generate(const Vector2f &sample, const hit_record &hrec) const override
    {
        Vector3f R = reflect(-hrec.wi, uvw.w());

        /* Sample from a Phong lobe centered around (0, 0, 1) */
        float sinAlpha = std::sqrt(1 - std::pow(sample[1], 2 / (specular_exponent + 1)));
        float cosAlpha = std::pow(sample[1], 1 / (specular_exponent + 1));
        float phi = (2.0f * M_PI) * sample[0];
        Vector3f localDir = Vector3f(
            sinAlpha * std::cos(phi),
            sinAlpha * std::sin(phi),
            cosAlpha
        );

        return onb(R).fromLocal(localDir);
    }

    onb uvw;
    const float specular_exponent;
    const ray r_in;
};

class constant_pdf : public pdf
{
public:
    constant_pdf(const float &p) : p(p) { }
    float value(const hit_record &hrec, const Vector3f &direction) const override
    {
        return p;
    }
    
    Vector3f generate(const Vector2f &sample, const hit_record &hrec) const override
    {
        throw std::runtime_error("Not implemented!");
    }

    const float p;
};

class mixture_pdf : public pdf
{
public:
    mixture_pdf(std::unique_ptr<pdf> p0, std::unique_ptr<pdf> p1, const float &prob_pdf, const float &sample) : prob_pdf(prob_pdf), rand_var(sample)
    {
        p[0] = std::move(p0);
        p[1] = std::move(p1);
    }

    float value(const hit_record &rec, const Vector3f &direction) const override
    {
        return prob_pdf*p[0]->value(rec, direction) + (1 - prob_pdf)*p[1]->value(rec, direction);
    }

    Vector3f generate(const Vector2f &sample, const hit_record &hrec) const override
    {
        if (rand_var < prob_pdf)
            return p[0]->generate(sample, hrec);

        return p[1]->generate(sample, hrec);
    }

    std::unique_ptr<pdf> p[2];
    const float prob_pdf;
    const float rand_var;
};

/* Implementation of rough conductor BSDF PDF functions taken from mitsuba rendeer */
class roughconductor_pdf : public pdf
{
public:
    roughconductor_pdf(const ray &r_in, const Vector3f &w, const float alphaU_, const float alphaV_, microfacet_distributions distribution_type_)
        : r_in(r_in), alphaU(alphaU_), alphaV(alphaV_), distribution_type(distribution_type_)  { uvw.build_from_w(w); }

    float value(const hit_record &hrec, const Vector3f &wo) const override
    {
        if (dot(hrec.normal, wo) <= 0
         || dot(hrec.normal, hrec.wi) <= 0)
            return 0.0f;


        /* Calculate the reflection half-vector */
        Vector3f H = unit_vector(wo+hrec.wi);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        float cosTheta = dot(hrec.normal, H);

        if (cosTheta <= 0)
            return 0.0f;

        float cosTheta2 = cosTheta * cosTheta;
        float beckmannExponent = ((H[0]*H[0]) / (alphaU * alphaU)
                + (H[1] * H[1]) / (alphaV * alphaV)) / cosTheta2;

        float eval_result;
        switch (distribution_type)
        {
        case microfacet_distributions::beckmann:
        {
            /* Beckmann distribution function for Gaussian random surfaces - [Walter 2005] evaluation */
            eval_result = fastexp(-beckmannExponent) /
                (M_PI * alphaU * alphaV * cosTheta2 * cosTheta2);
            break;
        }
        case microfacet_distributions::ggx:
        {
            /* GGX / Trowbridge-Reitz distribution function for rough surfaces */
            float root = ((float)1 + beckmannExponent) * cosTheta2;
            eval_result = (float)1 / (M_PI * alphaU * alphaV * root * root);
            break;
        }
        default:
            // invalid microfacet distribution type
            return -1;
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (eval_result * cosTheta < 1e-20f)
            eval_result = 0.0f;

        return eval_result * microfacet::smithG1(hrec.wi, H, hrec, alphaU, alphaV, distribution_type)
            / (4.0f * dot(hrec.wi, hrec.normal));
    }

    /* Taken from mitsuba */
    Vector3f fresnelConductorExact(float cosThetaI, const Vector3f &eta, const Vector3f &k) const
    {
        /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

        float cosThetaI2 = cosThetaI*cosThetaI,
              sinThetaI2 = 1-cosThetaI2,
              sinThetaI4 = sinThetaI2*sinThetaI2;

        Vector3f temp1 = eta*eta - k*k - Vector3f(sinThetaI2, sinThetaI2, sinThetaI2),
                 a2pb2 = safe_sqrt(temp1*temp1 + k*k*eta*eta*4),
                 a     = safe_sqrt((a2pb2 + temp1) * 0.5f);

        Vector3f term1 = a2pb2 + Vector3f(cosThetaI2, cosThetaI2, cosThetaI2),
                 term2 = a*(2*cosThetaI);

        Vector3f Rs2 = (term1 - term2) / (term1 + term2);

        Vector3f term3 = a2pb2*cosThetaI2 + Vector3f(sinThetaI4, sinThetaI4, sinThetaI4),
                 term4 = term2*sinThetaI2;

        Vector3f Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

        return 0.5f * (Rp2 + Rs2);
    }

    float pdfVisible(const Vector3f &wi, const Vector3f &m, const hit_record &hrec) const
    {
        float cosTheta = dot(hrec.normal, wi);
        if (cosTheta == 0)
            return 0.0f;
        return microfacet::smithG1(wi, m, hrec, alphaU, alphaV, distribution_type)
                * std::abs(dot(wi, m)) * microfacet::eval(m, hrec, alphaU, alphaV, distribution_type) / std::abs(cosTheta);
    }

    Vector2f sampleVisible11(float thetaI, Vector2f sample) const {
        const float SQRT_PI_INV = 1 / std::sqrt(M_PI);
        Vector2f slope;

        switch (distribution_type) {
        case  microfacet_distributions::beckmann: 
        {
            /* Special case (normal incidence) */
            if (thetaI < 1e-4f) {
                float sinPhi, cosPhi;
                float r = std::sqrt(-fastlog(1.0f - sample.x));
                sincos(2 * M_PI * sample.y, &sinPhi, &cosPhi);
                return Vector2(r * cosPhi, r * sinPhi);
            }

            /* The original inversion routine from the paper contained
               discontinuities, which causes issues for QMC integration
               and techniques like Kelemen-style MLT. The following code
               performs a numerical inversion with better behavior */
            float tanThetaI = std::tan(thetaI);
            float cotThetaI = 1 / tanThetaI;

            /* Search interval -- everything is parameterized
               in the erf_() domain */
            float a = -1, c = erf_(cotThetaI);
            float sample_x = std::max(sample.x, (float)1e-6f);

            /* We can do better (inverse of an approximation computed in Mathematica) */
            float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
            float b = c - (1 + c) * std::pow(1 - sample_x, fit);

            /* Normalization factor for the CDF */
            float normalization = 1 / (1 + c + SQRT_PI_INV *
                tanThetaI * fastexp(-cotThetaI * cotThetaI));

            int it = 0;
            while (++it < 10) {
                /* Bisection criterion -- the oddly-looking
                   boolean expression are intentional to check
                   for NaNs at little additional cost */
                if (!(b >= a && b <= c))
                    b = 0.5f * (a + c);

                /* Evaluate the CDF and its derivative
                   (i.e. the density function) */
                float invErf = erfinv(b);
                float value = normalization * (1 + b + SQRT_PI_INV *
                    tanThetaI * fastexp(-invErf * invErf)) - sample_x;
                float derivative = normalization * (1
                    - invErf * tanThetaI);

                if (std::abs(value) < 1e-5f)
                    break;

                /* Update bisection intervals */
                if (value > 0)
                    c = b;
                else
                    a = b;

                b -= value / derivative;
            }

            /* Now convert back into a slope value */
            slope.x = erfinv(b);

            /* Simulate Y component */
            slope.y = erfinv(2.0f * std::max(sample.y, (float)1e-6f) - 1.0f);
        };
        break;

        case  microfacet_distributions::ggx: 
        {
            /* Special case (normal incidence) */
            if (thetaI < 1e-4f) {
                float sinPhi, cosPhi;
                float r = safe_sqrt(sample.x / (1 - sample.x));
                sincos(2 * M_PI * sample.y, &sinPhi, &cosPhi);
                return Vector2(r * cosPhi, r * sinPhi);
            }

            /* Precomputations */
            float tanThetaI = std::tan(thetaI);
            float a = 1 / tanThetaI;
            float G1 = 2.0f / (1.0f + safe_sqrt(1.0f + 1.0f / (a * a)));

            /* Simulate X component */
            float A = 2.0f * sample.x / G1 - 1.0f;
            if (std::abs(A) == 1)
                A -= _copysign(1.0f, A) * EPSILON;
            float tmp = 1.0f / (A * A - 1.0f);
            float B = tanThetaI;
            float D = safe_sqrt(B * B * tmp * tmp - (A * A - B * B) * tmp);
            float slope_x_1 = B * tmp - D;
            float slope_x_2 = B * tmp + D;
            slope.x = (A < 0.0f || slope_x_2 > 1.0f / tanThetaI) ? slope_x_1 : slope_x_2;

            /* Simulate Y component */
            float S;
            if (sample.y > 0.5f) {
                S = 1.0f;
                sample.y = 2.0f * (sample.y - 0.5f);
            }
            else {
                S = -1.0f;
                sample.y = 2.0f * (0.5f - sample.y);
            }

            /* Improved fit */
            float z =
                (sample.y * (sample.y * (sample.y * (-(float)0.365728915865723) + (float)0.790235037209296) -
                (float)0.424965825137544) + (float)0.000152998850436920) /
                    (sample.y * (sample.y * (sample.y * (sample.y * (float)0.169507819808272 - (float)0.397203533833404) -
                (float)0.232500544458471) + (float)1) - (float)0.539825872510702);

            slope.y = S * z * std::sqrt(1.0f + slope.x * slope.x);
        };
                 break;

        default:
            return Vector2f(-1, -1);
        };
        return slope;
    }

    inline Vector3f sampleVisible(const Vector3f &_wi, const Vector2f &sample) const
    {
        /* Step 1: stretch wi */
        Vector3f wi = unit_vector(Vector3f(
            alphaU * _wi[0],
            alphaV * _wi[1],
            _wi[2]
        ));

        /* Get polar coordinates */
        float theta = 0, phi = 0;
        if (wi[2] < (float)0.99999) 
        {
            theta = std::acos(wi[2]);
            phi = std::atan2(wi[1], wi[0]);
        }
        float sinPhi, cosPhi;
        sincos(phi, &sinPhi, &cosPhi);

        /* Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1) */
        Vector2 slope = sampleVisible11(theta, sample);

        assert(std::isfinite(slope[0]) && std::isfinite(slope[1]));

        /* Step 3: rotate */
        slope = Vector2(
            cosPhi * slope.x - sinPhi * slope.y,
            sinPhi * slope.x + cosPhi * slope.y);

        /* Step 4: unstretch */
        slope.x *= alphaU;
        slope.y *= alphaV;

        /* Step 5: compute normal */
        float normalization = (float)1 / std::sqrt(slope.x * slope.x
            + slope.y * slope.y + (float)1.0);

        Vector3f result = Vector3f(
            -slope.x * normalization,
            -slope.y * normalization,
            normalization
        );

        assert(std::isfinite(result[0]) && std::isfinite(result[1]) && std::isfinite(result[2]));

        return result;
    }


    inline Vector3f sample_microfacet(const Vector3f& wi, const Vector2f& sample, float& pdf, const hit_record &hrec) const
    {
        Vector3f m;
        m = sampleVisible(wi, sample);
        //pdf = pdfVisible(wi, m, hrec);
        pdf = 0.0f;
        return m;
    }

    Vector3f generate(const Vector2f &sample, const hit_record &hrec) const override
    {

        //if (dot(hrec.normal, hrec.wi) < 0)
        //    return Vector3f(0.0f, 0.0f, 0.0f);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */

        /* Sample M, the microfacet normal */
        float pdf;
        Vector3f m = sample_microfacet(hrec.wi, sample, pdf, hrec);

        //if (pdf == 0)
        //    return Vector3f(0.0f, 0.0f, 0.0f);

        /* Perfect specular reflection based on the microfacet normal */
        Vector3f wo = reflect(hrec.wi, m);

        return wo;
    }

    onb uvw;
    const ray r_in;
    float alphaU, alphaV;
    microfacet_distributions distribution_type;
};

#endif
