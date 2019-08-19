#include "path_prt.h"
#include "prt.h"

using namespace PRT;

path_prt::path_prt(Scene *scene, const int &n_samples) : scene(scene), n_samples(n_samples)
{
	samples = PreComputeSamples(n_samples, n_bands);
	SH_project_environment();
}

// Here, n_coeffs = n_bands*n_bands and n_samples = sqrt_n_samples*sqrt_n_samples
void path_prt::SH_project_environment()
{
	Li_coeffs.reset(new Vector3f[n_coeffs]);
	std::memset(Li_coeffs.get(), 0, n_coeffs * sizeof(Vector3f));
	for (int i = 0; i < n_coeffs; ++i) {

		// For each sample
		for (int i = 0; i < n_samples; ++i) {
			double theta = samples[i].theta;
			double phi = samples[i].phi;
			for (int n = 0; n < n_coeffs; ++n) {
				Li_coeffs[n] += scene->env_map->eval(theta, phi) * samples[i].Ylm[n];
			}
		}
		// Divide the result by weight and number of samples
		const double factor = 4.0 * M_PI / n_samples;
		for (int i = 0; i < n_coeffs; ++i) {
			Li_coeffs[i] *= factor;
		}
	}
}

// TODO
// Convert from recursive to iterative
Vector3f path_prt::Li(const ray &r, Scene *scene, const int &depth, const hit_record &prev_hrec,
    const float &prev_bsdf_pdf, const int &s)
{
	hit_record hrec;
	auto &world = scene->world;
	auto &lights = scene->lights;

	if (world->hit(r, EPSILON, FLT_MAX, hrec))
	{
		scatter_record srec;
		if (hrec.mat_ptr->scatter(r, hrec, srec))
		{
			/* Sample BSDF to generate next ray direction for indirect lighting */
			hrec.p = hrec.p + (EPSILON * hrec.normal);
			//ray wo = srec.is_specular ? srec.specular_ray : ray(hrec.p, srec.pdf_ptr->generate());
			ray wo(hrec.p, unit_vector(samples[s].direction));
			const float surface_bsdf_pdf = srec.pdf_ptr->value(hrec, wo.direction());
			const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(wo, hrec, wo.direction());
			/* Reject current path in case the ray is on the wrong side of the surface (BRDF is 0 as ray is pointing away from the hemisphere )*/
			//if (surface_bsdf_pdf == 0)
			//{
			//	return Vector3f(0, 0, 0);
			//}
			const float cos_wi = std::max(dot(hrec.normal, wo.direction()), 0.0f);
			Vector3f Li(0.0f, 0.0f, 0.0f);
			if (cos_wi > 0.0f)
			{
				for (int i = 0; i < n_coeffs; ++i) {
					Li += Li_coeffs[i] * surface_bsdf * cos_wi * samples[s].Ylm[i];
				}
			}
			return Li;
			//return surface_bsdf * Li(wo, scene, depth + 1, hrec, surface_bsdf_pdf) * cos_wi / surface_bsdf_pdf;
		}
	}
	return scene->env_map->eval(r, hrec, depth);
}

