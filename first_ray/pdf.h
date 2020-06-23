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
    const double &r0 = sample[0];
    const double &r1 = sample[1];
    const double r = sqrt(r0);
    const double phi = 2 * (double)M_PI * r1;
    const double x = r * cos(phi);
    const double y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}

inline Vector3f hemisphere_to_cosine_direction(double &theta, double& phi, const Vector2f &sample)
{
    const double &r0 = sample[0];
    const double &r1 = sample[1];
    const double r = sqrt(r0);
    theta = asin(r);
    phi = 2 * (double)M_PI * r1;
    const double x = r * cos(phi);
    const double y = r * sin(phi);

    return Vector3f(x, y, sqrt(1 - r0));
}

inline Vector3f uniform_sample_sphere(const Vector2f &sample) {
    Vector2f u(sample);
    const double  z = 1 - 2 * u[0];
    const double r = std::sqrt(std::max((double)0, (double)1 - z * z));
    const double phi = 2 * (double)M_PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline Vector3f random_to_sphere(const double &radius, const double &distance_squared, const Vector2f& sample)
{
    const double &r1 = sample[0];
    const double &r2 = sample[1];
    double z = 1 + r2 * (sqrt(1 - radius*radius/distance_squared) - 1);
    double phi = 2*(double)M_PI*r1;
    double x = cos(phi)*sqrt(1-z*z);
    double y = sin(phi)*sqrt(1-z*z);

    return Vector3f(x, y, z);
}

class pdf;

struct scatter_record
{
    ray specular_ray;
    bool is_specular;
    hit_record hrec;
    double eta;
    std::unique_ptr<pdf> pdf_ptr;

    scatter_record(const hit_record &hrec_) : hrec(hrec_) {}
};

class pdf
{
public:
    virtual double value(const hit_record &, const Vector3f &) const = 0;
    virtual Vector3f generate(const Vector2f &sample, scatter_record &srec) const = 0;
    virtual ~pdf() = 0;
};

class cosine_pdf : public pdf
{
public:
    cosine_pdf(const Vector3f &w) { uvw.build_from_w(w); }
    double value(const hit_record &hrec, const Vector3f &direction) const override
    {
        double cosine = std::max<double>(dot(uvw.w(), unit_vector(direction)), 0.0);
        
        return cosine / (double) (double)M_PI;
    }

    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        return uvw.fromLocal(hemisphere_to_cosine_direction(sample));
    }

    onb uvw;
};

class cosine_power_pdf : public pdf
{
public:
    cosine_power_pdf(const ray &r_in, const Vector3f &w, const double &specular_exponent) : specular_exponent(specular_exponent), r_in(r_in) { uvw.build_from_w(w); }
    // TODO: modify this stuff - modified phong BRDF paper
    double value(const hit_record& hrec, const Vector3f& wo) const override
    {
        if (dot(hrec.normal, wo) <= 0
            || dot(hrec.normal, hrec.wi) <= 0)
            return 0.0;
        const double alpha = std::max(0.0, dot(reflect(-hrec.wi, uvw.w()), wo));
        const double specular_pdf = std::pow(alpha, (double)specular_exponent);

        return specular_pdf * (specular_exponent + 1.0) / (2 * (double)M_PI);
    }

    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        hit_record &hrec = srec.hrec;
        Vector3f R = reflect(-hrec.wi, uvw.w());

        /* Sample from a Phong lobe centered around (0, 0, 1) */
        double sinAlpha = std::sqrt(1 - std::pow(sample[1], 2 / (specular_exponent + 1)));
        double cosAlpha = std::pow(sample[1], 1 / (specular_exponent + 1));
        double phi = (2.0 * M_PI) * sample[0];
        Vector3f localDir = Vector3f(
            sinAlpha * std::cos(phi),
            sinAlpha * std::sin(phi),
            cosAlpha
        );

        return onb(R).fromLocal(localDir);
    }

    onb uvw;
    const double specular_exponent;
    const ray r_in;
};

class dielectric_pdf : public pdf
{
public:
    dielectric_pdf(const Vector3f &w, const double ref_idx_) : ref_idx(ref_idx_) { uvw.build_from_w(w); }
    double value(const hit_record &hrec, const Vector3f &wo) const override
    {
        double cosThetaT;
        //double cosWi = dot(hrec.wi, hrec.normal);
        double F = fresnelDielectricExt(dot(hrec.wi, hrec.normal), cosThetaT, ref_idx);

        if (dot(hrec.wi, hrec.normal) * dot(wo, hrec.normal) >= 0)
        {
            if (std::abs(dot(reflect(-hrec.wi, hrec.normal), wo) - 1) > DELTA_EPSILON)
                return 0.0;

            return F;
        }
        else
        {
            if (std::abs(dot(refract(hrec.wi, hrec.normal, ref_idx, cosThetaT), wo)-1) > DELTA_EPSILON)
                return 0.0;

            return 1.0 - F;
        }
    }

    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        hit_record &hrec = srec.hrec;
        double cosThetaT;
        double F = fresnelDielectricExt(dot(hrec.wi, hrec.normal), cosThetaT, ref_idx);

        if (sample.x <= F)
        {
            srec.eta = 1.0;
            return reflect(-hrec.wi, hrec.normal);
        }
        else
        {
            srec.eta = cosThetaT < 0 ? ref_idx : (1.0 / ref_idx);
            return refract(hrec.wi, hrec.normal, ref_idx, cosThetaT);
        }
    }

    double ref_idx;
    onb uvw;
};

class constant_pdf : public pdf
{
public:
    constant_pdf(const double &p) : p(p) { }
    double value(const hit_record &hrec, const Vector3f &direction) const override
    {
        return p;
    }
    
    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        throw std::runtime_error("Not implemented!");
    }

    const double p;
};

class mixture_pdf : public pdf
{
public:
    mixture_pdf(std::unique_ptr<pdf> p0, std::unique_ptr<pdf> p1, const double &prob_pdf, const double &sample) : prob_pdf(prob_pdf), rand_var(sample)
    {
        p[0] = std::move(p0);
        p[1] = std::move(p1);
    }

    double value(const hit_record &rec, const Vector3f &direction) const override
    {
        return prob_pdf*p[0]->value(rec, direction) + (1 - prob_pdf)*p[1]->value(rec, direction);
    }

    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        if (rand_var < prob_pdf)
            return p[0]->generate(sample, srec);

        return p[1]->generate(sample, srec);
    }

    std::unique_ptr<pdf> p[2];
    const double prob_pdf;
    const double rand_var;
};

/* Implementation of rough conductor BSDF PDF functions taken from mitsuba rendeer */
class roughconductor_pdf : public pdf
{
public:
    roughconductor_pdf(const ray &r_in, const Vector3f &w, const double alphaU_, const double alphaV_, microfacet_distributions distribution_type_)
        : r_in(r_in), alphaU(alphaU_), alphaV(alphaV_), distribution_type(distribution_type_)  { uvw.build_from_w(w); }

    double value(const hit_record &hrec, const Vector3f &wo) const override
    {
        if (dot(hrec.normal, wo) <= 0
         || dot(hrec.normal, hrec.wi) <= 0)
            return 0.0;


        /* Calculate the reflection half-vector */
        Vector3f H = unit_vector(wo+hrec.wi);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        double cosTheta = dot(hrec.normal, H);

        if (cosTheta <= 0)
            return 0.0;

        double cosTheta2 = cosTheta * cosTheta;
        double beckmannExponent = ((H[0]*H[0]) / (alphaU * alphaU)
                + (H[1] * H[1]) / (alphaV * alphaV)) / cosTheta2;

        double eval_result;
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
                double root = ((double)1 + beckmannExponent) * cosTheta2;
                eval_result = (double)1 / (M_PI * alphaU * alphaV * root * root);
                break;
            }
            default:
                // invalid microfacet distribution type
                return -1;
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (eval_result * cosTheta < 1e-20f)
            eval_result = 0.0;

        return eval_result * microfacet::smithG1(hrec.wi, H, hrec, alphaU, alphaV, distribution_type)
            / (4.0 * dot(hrec.wi, hrec.normal));
    }

    /* Taken from mitsuba */
    Vector3f fresnelConductorExact(double cosThetaI, const Vector3f &eta, const Vector3f &k) const
    {
        /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

        double cosThetaI2 = cosThetaI*cosThetaI,
              sinThetaI2 = 1-cosThetaI2,
              sinThetaI4 = sinThetaI2*sinThetaI2;

        Vector3f temp1 = eta*eta - k*k - Vector3f(sinThetaI2, sinThetaI2, sinThetaI2),
                 a2pb2 = safe_sqrt(temp1*temp1 + k*k*eta*eta*4),
                 a     = safe_sqrt((a2pb2 + temp1) * 0.5);

        Vector3f term1 = a2pb2 + Vector3f(cosThetaI2, cosThetaI2, cosThetaI2),
                 term2 = a*(2*cosThetaI);

        Vector3f Rs2 = (term1 - term2) / (term1 + term2);

        Vector3f term3 = a2pb2*cosThetaI2 + Vector3f(sinThetaI4, sinThetaI4, sinThetaI4),
                 term4 = term2*sinThetaI2;

        Vector3f Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

        return 0.5 * (Rp2 + Rs2);
    }

    double pdfVisible(const Vector3f &wi, const Vector3f &m, const hit_record &hrec) const
    {
        double cosTheta = dot(hrec.normal, wi);
        if (cosTheta == 0)
            return 0.0;
        return microfacet::smithG1(wi, m, hrec, alphaU, alphaV, distribution_type)
                * std::abs(dot(wi, m)) * microfacet::eval(m, hrec, alphaU, alphaV, distribution_type) / std::abs(cosTheta);
    }

    Vector2f sampleVisible11(double thetaI, Vector2f sample) const {
        const double SQRT_PI_INV = 1 / std::sqrt(M_PI);
        Vector2f slope;

        switch (distribution_type) {
        case  microfacet_distributions::beckmann: 
        {
            /* Special case (normal incidence) */
            if (thetaI < 1e-4f) {
                double sinPhi, cosPhi;
                double r = std::sqrt(-fastlog(1.0 - sample.x));
                sincos(2 * M_PI * sample.y, &sinPhi, &cosPhi);
                return Vector2(r * cosPhi, r * sinPhi);
            }

            /* The original inversion routine from the paper contained
               discontinuities, which causes issues for QMC integration
               and techniques like Kelemen-style MLT. The following code
               performs a numerical inversion with better behavior */
            double tanThetaI = std::tan(thetaI);
            double cotThetaI = 1 / tanThetaI;

            /* Search interval -- everything is parameterized
               in the erf_() domain */
            double a = -1, c = erf_(cotThetaI);
            double sample_x = std::max(sample.x, (double)1e-6f);

            /* We can do better (inverse of an approximation computed in Mathematica) */
            double fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
            double b = c - (1 + c) * std::pow(1 - sample_x, fit);

            /* Normalization factor for the CDF */
            double normalization = 1 / (1 + c + SQRT_PI_INV *
                tanThetaI * std::exp(-cotThetaI * cotThetaI));

            int it = 0;
            while (++it < 10) {
                /* Bisection criterion -- the oddly-looking
                   boolean expression are intentional to check
                   for NaNs at little additional cost */
                if (!(b >= a && b <= c))
                    b = 0.5 * (a + c);

                /* Evaluate the CDF and its derivative
                   (i.e. the density function) */
                double invErf = erfinv(b);
                double value = normalization * (1 + b + SQRT_PI_INV *
                    tanThetaI * fastexp(-invErf * invErf)) - sample_x;
                double derivative = normalization * (1
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
            slope.y = erfinv(2.0 * std::max(sample.y, (double)1e-6f) - 1.0);
        };
        break;

        case microfacet_distributions::ggx:
        {
            /* Special case (normal incidence) */
            if (thetaI < 1e-4f) {
                double sinPhi, cosPhi;
                double r = safe_sqrt(sample.x / (1 - sample.x));
                sincos(2 * M_PI * sample.y, &sinPhi, &cosPhi);
                return Vector2(r * cosPhi, r * sinPhi);
            }

            /* Precomputations */
            double tanThetaI = std::tan(thetaI);
            double a = 1 / tanThetaI;
            double G1 = 2.0 / (1.0 + safe_sqrt(1.0 + 1.0 / (a * a)));

            /* Simulate X component */
            double A = 2.0 * sample.x / G1 - 1.0;
            if (std::abs(A) == 1)
                A -= copysignf(1.0, A) * EPSILON;
            double tmp = 1.0 / (A * A - 1.0);
            double B = tanThetaI;
            double D = safe_sqrt(B * B * tmp * tmp - (A * A - B * B) * tmp);
            double slope_x_1 = B * tmp - D;
            double slope_x_2 = B * tmp + D;
            slope.x = (A < 0.0 || slope_x_2 > 1.0 / tanThetaI) ? slope_x_1 : slope_x_2;

            /* Simulate Y component */
            double S;
            if (sample.y > 0.5) {
                S = 1.0;
                sample.y = 2.0 * (sample.y - 0.5);
            }
            else {
                S = -1.0;
                sample.y = 2.0 * (0.5 - sample.y);
            }

            /* Improved fit */
            double z =
                (sample.y * (sample.y * (sample.y * (-(double)0.365728915865723) + (double)0.790235037209296) -
                (double)0.424965825137544) + (double)0.000152998850436920) /
                    (sample.y * (sample.y * (sample.y * (sample.y * (double)0.169507819808272 - (double)0.397203533833404) -
                (double)0.232500544458471) + (double)1) - (double)0.539825872510702);

            slope.y = S * z * std::sqrt(1.0 + slope.x * slope.x);
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
        double theta = 0, phi = 0;
        if (wi[2] < (double)0.99999) 
        {
            theta = std::acos(wi[2]);
            phi = std::atan2(wi[1], wi[0]);
        }
        double sinPhi, cosPhi;
        sincos(phi, &sinPhi, &cosPhi);

        /* Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1) */
        Vector2f slope = sampleVisible11(theta, sample);

        //assert(std::isfinite(slope[0]) && std::isfinite(slope[1]));
        if (!std::isfinite(slope.x))
            slope[0] = 0.0;

        /* Step 3: rotate */
        slope = Vector2(
            cosPhi * slope.x - sinPhi * slope.y,
            sinPhi * slope.x + cosPhi * slope.y);

        /* Step 4: unstretch */
        slope.x *= alphaU;
        slope.y *= alphaV;

        /* Step 5: compute normal */
        double normalization = (double)1 / std::sqrt(slope.x * slope.x
            + slope.y * slope.y + (double)1.0);

        Vector3f result = Vector3f(
            -slope.x * normalization,
            -slope.y * normalization,
            normalization
        );

        assert(std::isfinite(result[0]) && std::isfinite(result[1]) && std::isfinite(result[2]));

        return result;
    }

    inline Vector3f sample_microfacet(const Vector3f& wi, const Vector2f& sample, double& pdf, const hit_record &hrec) const
    {
        Vector3f m;
        m = sampleVisible(wi, sample);
        //pdf = pdfVisible(wi, m, hrec);
        pdf = 0.0;
        return m;
    }

    Vector3f generate(const Vector2f &sample, scatter_record &srec) const override
    {
        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */

        hit_record &hrec = srec.hrec;
        /* Sample M, the microfacet normal */
        double pdf;
        Vector3f m = sample_microfacet(uvw.toLocal(hrec.wi), sample, pdf, hrec);

        /* Perfect specular reflection based on the microfacet normal */
        Vector3f wo = reflect(-hrec.wi, uvw.fromLocal(m));

        return wo;
    }

    onb uvw;
    const ray r_in;
    double alphaU, alphaV;
    microfacet_distributions distribution_type;
};

#endif
