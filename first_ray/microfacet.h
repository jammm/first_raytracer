#ifndef MICROFACET_H
#define MICROFACET_H

#include "geometry.h"

struct microfacet
{
    static Vector3f fresnelConductorExact(double cosThetaI, const Vector3f &eta, const Vector3f &k)
    {
        /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

        double cosThetaI2 = cosThetaI * cosThetaI,
            sinThetaI2 = 1 - cosThetaI2,
            sinThetaI4 = sinThetaI2 * sinThetaI2;

        Vector3f temp1 = eta * eta - k * k - Vector3f(sinThetaI2, sinThetaI2, sinThetaI2),
            a2pb2 = safe_sqrt(temp1 * temp1 + k * k * eta * eta * 4),
            a = safe_sqrt((a2pb2 + temp1) * 0.5);

        Vector3f term1 = a2pb2 + Vector3f(cosThetaI2, cosThetaI2, cosThetaI2),
            term2 = a * (2 * cosThetaI);

        Vector3f Rs2 = (term1 - term2) / (term1 + term2);

        Vector3f term3 = a2pb2 * cosThetaI2 + Vector3f(sinThetaI4, sinThetaI4, sinThetaI4),
            term4 = term2 * sinThetaI2;

        Vector3f Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

        return 0.5 * (Rp2 + Rs2);
    }

    static double projectRoughness(const Vector3f& v, const hit_record& hrec, const double alphaU, const double alphaV)
    {
        double cosTheta = dot(hrec.normal, v);
        double invSinTheta2 = 1 / (1 - cosTheta * cosTheta);

        if ((alphaU == alphaV) || invSinTheta2 <= 0)
            return alphaU;

        double cosPhi2 = v[0] * v[0] * invSinTheta2;
        double sinPhi2 = v[1] * v[1] * invSinTheta2;

        return std::sqrt(cosPhi2 * alphaU * alphaU + sinPhi2 * alphaV * alphaV);
    }

	static double smithG1(const Vector3f &v, const Vector3f &m, const hit_record &hrec, const double alphaU, const double alphaV, microfacet_distributions distribution_type)
    {
        const double cosTheta = dot(hrec.normal, v);
        /* Ensure consistent orientation (can't see the back
            of the microfacet from the front and vice versa) */
        if (dot(v, m) * cosTheta <= 0)
            return 0.0;

        /* Perpendicular incidence -- no shadowing/masking */
        double temp = 1 - (cosTheta * cosTheta);
        if (temp <= 0.0)
            return 1.0;
        const double tanTheta = std::sqrt(temp) / cosTheta;

        double alpha = projectRoughness(v, hrec, alphaU, alphaV);
        switch (distribution_type)
        {
            case microfacet_distributions::beckmann:
            {
                double a = 1.0 / (alpha * tanTheta);
                if (a >= 1.6f)
                    return 1.0;

                /* Use a fast and accurate (<0.35% rel. error) rational
                    approximation to the shadowing-masking function */
                double aSqr = a * a;
                return (3.535f * a + 2.181f * aSqr)
                    / (1.0 + 2.276f * a + 2.577f * aSqr);
                break;
            }

            case microfacet_distributions::ggx:
            {
                double root = alpha * tanTheta;
                return 2.0 / (1.0 + hypot2(1.0, root));
                break;
            }

            default:
                return -1.0;
        }
    }

    static double eval(const Vector3f &m, const hit_record &hrec, const double alphaU, const double alphaV, microfacet_distributions distribution_type)
    {
        Vector3f m_ = onb(hrec.normal).toLocal(m);
        const double cosTheta = m_.z();
        if (cosTheta <= 0)
            return 0.0;

        double cosTheta2 = cosTheta * cosTheta;
        double beckmannExponent = ((m_[0] * m_[0]) / (alphaU * alphaU)
            + (m_[1] * m_[1]) / (alphaV * alphaV)) / cosTheta2;

        double result;
        switch (distribution_type) {
        case microfacet_distributions::beckmann:
        {
            /* Beckmann distribution function for Gaussian random surfaces - [Walter 2005] evaluation */
            result = fastexp(-beckmannExponent) /
                (M_PI * alphaU * alphaV * cosTheta2 * cosTheta2);
            break;
        }

        case microfacet_distributions::ggx:
        {
            /* GGX / Trowbridge-Reitz distribution function for rough surfaces */
            double root = ((double)1 + beckmannExponent) * cosTheta2;
            result = (double)1 / (M_PI * alphaU * alphaV * root * root);
            break;
        }

        default:
            return -1;
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (result * cosTheta < 1e-20f)
            result = 0;

        return result;
    }
};


#endif