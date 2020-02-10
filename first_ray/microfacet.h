#ifndef MICROFACET_H
#define MICROFACET_H

#include "geometry.h"

struct microfacet
{
    static Vector3f fresnelConductorExact(float cosThetaI, const Vector3f &eta, const Vector3f &k)
    {
        /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

        float cosThetaI2 = cosThetaI * cosThetaI,
            sinThetaI2 = 1 - cosThetaI2,
            sinThetaI4 = sinThetaI2 * sinThetaI2;

        Vector3f temp1 = eta * eta - k * k - Vector3f(sinThetaI2, sinThetaI2, sinThetaI2),
            a2pb2 = safe_sqrt(temp1 * temp1 + k * k * eta * eta * 4),
            a = safe_sqrt((a2pb2 + temp1) * 0.5f);

        Vector3f term1 = a2pb2 + Vector3f(cosThetaI2, cosThetaI2, cosThetaI2),
            term2 = a * (2 * cosThetaI);

        Vector3f Rs2 = (term1 - term2) / (term1 + term2);

        Vector3f term3 = a2pb2 * cosThetaI2 + Vector3f(sinThetaI4, sinThetaI4, sinThetaI4),
            term4 = term2 * sinThetaI2;

        Vector3f Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

        return 0.5f * (Rp2 + Rs2);
    }

    static float projectRoughness(const Vector3f& v, const hit_record& hrec, const float alphaU, const float alphaV)
    {
        float cosTheta = dot(hrec.normal, v);
        float invSinTheta2 = 1 / (1 - cosTheta * cosTheta);

        if ((alphaU == alphaV) || invSinTheta2 <= 0)
            return alphaU;

        float cosPhi2 = v[0] * v[0] * invSinTheta2;
        float sinPhi2 = v[1] * v[1] * invSinTheta2;

        return std::sqrt(cosPhi2 * alphaU * alphaU + sinPhi2 * alphaV * alphaV);
    }

	static float smithG1(const Vector3f &v, const Vector3f &m, const hit_record &hrec, const float alphaU, const float alphaV, microfacet_distributions distribution_type)
    {
        const float cosTheta = dot(hrec.normal, v);
        /* Ensure consistent orientation (can't see the back
            of the microfacet from the front and vice versa) */
        if (dot(v, m) * cosTheta <= 0)
            return 0.0f;

        /* Perpendicular incidence -- no shadowing/masking */
        float temp = 1 - (cosTheta * cosTheta);
        if (temp <= 0.0f)
            return 1.0f;
        const float tanTheta = std::sqrt(temp) / cosTheta;

        float alpha = projectRoughness(v, hrec, alphaU, alphaV);
        switch (distribution_type)
        {
            case microfacet_distributions::beckmann:
            {
                float a = 1.0f / (alpha * tanTheta);
                if (a >= 1.6f)
                    return 1.0f;

                /* Use a fast and accurate (<0.35% rel. error) rational
                    approximation to the shadowing-masking function */
                float aSqr = a * a;
                return (3.535f * a + 2.181f * aSqr)
                    / (1.0f + 2.276f * a + 2.577f * aSqr);
                break;
            }

            case microfacet_distributions::ggx:
            {
                float root = alpha * tanTheta;
                return 2.0f / (1.0f + hypot2(1.0f, root));
                break;
            }

            default:
                return -1.0f;
        }
    }

    static float eval(const Vector3f &m, const hit_record &hrec, const float alphaU, const float alphaV, microfacet_distributions distribution_type)
    {
        const float cosTheta = dot(hrec.normal, m);
        if (cosTheta <= 0)
            return 0.0f;

        float cosTheta2 = cosTheta * cosTheta;
        float beckmannExponent = ((m[0] * m[0]) / (alphaU * alphaU)
            + (m[1] * m[1]) / (alphaV * alphaV)) / cosTheta2;

        float result;
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
            float root = ((float)1 + beckmannExponent) * cosTheta2;
            result = (float)1 / (M_PI * alphaU * alphaV * root * root);
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