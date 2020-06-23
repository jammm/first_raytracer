#ifndef MATERIAL_H
#define MATERIAL_H

#include "hitable.h"
#include "util.h"
#include "texture.h"
#include "onb.h"
#include "microfacet.h"
#include "pdf.h"
#include "image.h"

#include <algorithm>
#include <cctype>
#include <string>

inline Vector3f random_in_unit_sphere(const Vector3f &sample)
{
    Vector3f p;

    do
    {
        p = 2.0 * sample - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0);

    return p;
}

inline Vector3f random_on_unit_sphere(const Vector3f &sample)
{
    Vector3f p;

    do
    {
        p = 2.0 * sample - Vector3f(1, 1, 1);
    } while (p.squared_length() >= 1.0);

    return unit_vector(p);
}

class material
{
public:
    virtual bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const { return false; }
    virtual Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const { return Vector3f(0, 0, 0); }
    virtual Vector3f get_albedo(const hit_record& rec) const { return Vector3f(0, 0, 0); }
    virtual Vector3f emitted(const ray &r_in, const hit_record &rec) const { return Vector3f(0, 0, 0); }
    virtual ~material();
};

class lambertian : public material
{
public:
    lambertian(texture *a) : albedo(a) {}

    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override
    {
        srec.is_specular = false;
        srec.pdf_ptr = std::make_unique<cosine_pdf>(hrec.normal);
        return true;
    }

    Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const override
    {
        return albedo->value(rec) / M_PI;
    }

    Vector3f get_albedo(const hit_record& rec) const override
    {
        return albedo->value(rec);
    }

    texture *albedo;
};

class modified_phong : public material
{
public:
    modified_phong(texture *diffuse_reflectance_, texture *specular_reflectance_, const double &specular_exponent) 
        : diffuse_reflectance(diffuse_reflectance_), specular_reflectance(specular_reflectance_), specular_exponent(specular_exponent) {}

    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override
    {
        //const double rand_var = gen_cano_rand();
        //constexpr double specular_chance = 0.5;

        srec.is_specular = true;

        srec.pdf_ptr = std::make_unique<cosine_power_pdf>(r_in, hrec.normal, specular_exponent);
        srec.specular_ray = ray(hrec.p, srec.pdf_ptr->generate(Vector2f(sample[0], sample[1]), srec));

        return true;
    }

    Vector3f eval_bsdf(const ray &r_in, const hit_record &hrec, const Vector3f &wo) const override
    {
        const Vector3f reflected = reflect(-hrec.wi, hrec.normal);
        const double alpha = std::max(0.0, dot(reflected, wo));
        const Vector3f result = (diffuse_reflectance->value(hrec) / M_PI) + specular_reflectance->value(hrec)
                * ((specular_exponent + 2) / (2 * M_PI)) * std::pow(alpha, specular_exponent);

        return result * dot(hrec.normal, wo);
    }

    texture *diffuse_reflectance;
    texture *specular_reflectance;
    const double specular_exponent;
    Vector3f wo;
};

class metal : public material
{
public:
    metal(const Vector3f &a, double f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1.0; }
    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override
    {
        Vector3f reflected = reflect(unit_vector(r_in.direction()), hrec.normal);
        srec.specular_ray = ray(hrec.p, reflected + fuzz*random_in_unit_sphere(sample));
        srec.is_specular = true;
        srec.pdf_ptr = std::make_unique<constant_pdf>(1.0);
        return true;
    }
    Vector3f eval_bsdf(const ray &r_in, const hit_record &rec, const Vector3f &wo) const override
    {
        return albedo;
    }

    Vector3f albedo;
    double fuzz;
};

/* Dielectric BSDF taken from mitsuba */
class dielectric : public material
{
public:
    dielectric(double ref_idx_, texture *specular_reflectance_, texture *specular_transmittance_)
        : ref_idx(ref_idx_), specular_reflectance(specular_reflectance_), specular_transmittance(specular_transmittance_) {}

    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override
    {
        srec.is_specular = true;
        srec.pdf_ptr = std::make_unique<dielectric_pdf>(hrec.normal, ref_idx);
        srec.specular_ray = ray(hrec.p, srec.pdf_ptr->generate(Vector2f(sample[0], sample[1]), srec));

        return true;
    }
    Vector3f eval_bsdf(const ray &r_in, const hit_record &hrec, const Vector3f &wo) const override
    {
        //const double cosWi = dot(hrec.wi, hrec.normal);
        double cosThetaT;
        double F = fresnelDielectricExt(dot(hrec.wi, hrec.normal), cosThetaT, ref_idx);

        if (dot(hrec.wi, hrec.normal) * dot(wo, hrec.normal) >= 0)
        {
            if (std::abs(dot(reflect(-hrec.wi, hrec.normal), wo)-1) > DELTA_EPSILON)
                return Vector3f(0.0, 0.0, 0.0);

            return specular_reflectance->value(hrec) * F;
        } 
        else 
        {
            //onb uvw(hrec.normal);
            if (std::abs(dot(refract(hrec.wi, hrec.normal, ref_idx, cosThetaT), wo)-1) > DELTA_EPSILON)
                return Vector3f(0.0, 0.0, 0.0);

            /* Radiance must be scaled to account for the solid angle compression
               that occurs when crossing the interface. */
            double factor = cosThetaT < 0 ? (1.0/ref_idx) : (ref_idx);

            return specular_reflectance->value(hrec) * factor * factor * (1 - F);
        }
    }

    double ref_idx;
    texture *specular_reflectance;
    texture *specular_transmittance;
};

class diffuse_light : public material
{
public:
    diffuse_light(texture *a) : emit(a) {}
    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override { return false; }
    Vector3f emitted(const ray &r_in, const hit_record &rec) const override
    {
        if (dot(rec.normal, r_in.d) < 0)
            return emit->value(rec);
        else
            return Vector3f(0, 0, 0);
    }
    texture *emit;
};

class point_light_mat : public material
{
public:
    point_light_mat(texture* a) : emit(a) {}
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, const Vector3f& sample) const override { return false; }
    Vector3f emitted(const ray& r_in, const hit_record& rec) const override
    {
        return emit->value(rec) / (rec.p - r_in.o).squared_length();
    }
    texture* emit;
};

class environment_map : public material
{
public:
    environment_map(std::string env_map_filename) : env_map_filename(env_map_filename)
    {
        auto env_map_img = std::make_unique<image>(env_map_filename, formats::STBI_HDR);
        env_map_tex = std::make_unique<image_texture>(std::move(env_map_img));
    }
    environment_map(std::unique_ptr<texture> e)
    {
        env_map_tex = std::move(e);
    }
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, const Vector3f& sample) const override { return false; }
    Vector3f eval(const ray& r_in, hit_record hrec, const int &depth) const
    {
        const Vector3f direction = unit_vector(r_in.d);
        double phi = std::atan2(direction.x(), -direction.z());
        double theta = acos(direction.y());

        phi = (phi < 0) ? (phi + M_PI * 2) : phi;
        theta = (theta < 0) ? (theta + M_PI) : theta;

        hrec.u = phi / (2.0 * M_PI);
        hrec.v = theta / (M_PI);

        return env_map_tex->value(hrec);
    }
    Vector3f eval(const double &theta, const double &phi) const
    {
        hit_record rec;
        rec.u = (phi + M_PI*2) / (2.0*M_PI);
        rec.v = (theta + M_PI) / M_PI;

        return env_map_tex->value(rec);
    }

    std::string env_map_filename;
    std::unique_ptr<texture> env_map_tex;
};

class rough_conductor : public material
{
public:
    rough_conductor(const double alpha_, const double extEta_, const Vector3f &eta_, const Vector3f &k_,
                    texture *specular_reflectance_, const std::string &distribution_) 
        : alpha(alpha_), extEta(extEta_), alphaU(alpha), alphaV(alpha), eta(eta_), k(k_), specular_reflectance(specular_reflectance_), distribution(distribution_)
    {
        std::transform(distribution_.begin(), distribution_.end(), distribution.begin(), ::tolower);
        if (distribution == "ggx")
        {
            distribution_type = microfacet_distributions::ggx;
        }
        else
            distribution_type = microfacet_distributions::beckmann;
    }

    bool scatter(const ray &r_in, const hit_record &hrec, scatter_record &srec, const Vector3f &sample) const override
    {
        srec.is_specular = true;

        srec.pdf_ptr = std::make_unique<roughconductor_pdf>(r_in, hrec.normal, alphaU, alphaV, distribution_type);
        srec.specular_ray = ray(hrec.p, unit_vector(srec.pdf_ptr->generate(Vector2f(sample[0], sample[1]), srec)));

        return true;
    }

    inline double G_term(const Vector3f &wi, const Vector3f &wo, const Vector3f &m, const hit_record &hrec) const {
        return microfacet::smithG1(wi, m, hrec, alphaU, alphaV, distribution_type) * microfacet::smithG1(wo, m, hrec, alphaU, alphaV, distribution_type);
    }

    Vector3f eval_bsdf(const ray &r_in, const hit_record &hrec, const Vector3f &wo) const override
    {
        double cosWi = dot(hrec.wi, hrec.normal);
        /* Stop if this component was not requested */
        if (cosWi <= 0 ||
            dot(wo, hrec.normal) <= 0)
            return Vector3f(0.0, 0.0, 0.0);

        /* Calculate the reflection half-vector */
        Vector3f H = unit_vector(wo+hrec.wi);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */

        /* Evaluate the microfacet normal distribution */
        const double D = microfacet::eval(H, hrec, alphaU, alphaV, distribution_type);
        if (D == 0)
            return Vector3f(0.0, 0.0, 0.0);

        /* Fresnel factor */
        const Vector3f F = microfacet::fresnelConductorExact(dot(hrec.wi, H), eta, k) *
            specular_reflectance->value(hrec);

        /* Smith's shadow-masking function */
        const double G = G_term(hrec.wi, wo, H, hrec);

        /* Calculate the total amount of reflection */
        double model = D * G / (4.0 * cosWi);

        return F * model;
    }

    double alpha;
    double extEta;
    double alphaU, alphaV;
    Vector3f eta, k;
    texture *specular_reflectance;
    std::string distribution;
    microfacet_distributions distribution_type;
};

#endif
