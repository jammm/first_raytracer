#include "pssmlt.h"
#include "material.h"
#include <mutex>

// primary space Markov chain
inline double perturb(const float value, const float s1, const float s2, sampler &s) {
	double Result;
    double r = s.get1d();
	if (r < 0.5) {
		r = r * 2.0;
		Result = value + s2 * std::exp(-std::log(s2 / s1) * r); if (Result > 1.0f) Result -= 1.0f;
	} else {
		r = (r - 0.5) * 2.0;
		Result = value - s2 * std::exp(-std::log(s2 / s1) * r); if (Result < 0.0f) Result += 1.0f;
	}
	return Result;
}

void AccumulatePathContribution(const pssmlt::PathContribution pc, const float mScaling, viewer *v) {
    if (pc.sc == 0) return;
    const int ix = int(pc.x);
    const int iy = int(pc.y);
    const Vector3f c = pc.c * mScaling;
    const int PixelWidth = v->nx;
    const int PixelHeight = v->ny;
    const float scale = PixelWidth * PixelHeight / float(v->ns);
    static std::mutex m;

    if ((ix < 0) || (ix >= PixelWidth) || (iy < 0) || (iy >= PixelHeight))
        return;
    std::lock_guard<std::mutex> guard(m);
    for (int comp=0;comp<3;++comp)
    {
        v->fout_image[((ix + iy * PixelWidth) * 3) + comp] += scale * c[comp];
        v->out_image[((ix + iy * PixelWidth) * 3) + comp] = std::min(int(pow(v->fout_image[((ix + iy * PixelWidth) * 3) + comp] + c[comp] * scale, 1.0 / 2.2) * 255.9999), 255);
    }

}

struct TMarkovChain
{   
    float u[NumStates];
    pssmlt::PathContribution C;

    TMarkovChain()
    {
        sampler s;
        for (int i = 0; i < NumStates; i++) u[i] = s.get1d();
    }

	TMarkovChain(sampler &s)
    {
        for (int i = 0; i < NumStates; i++) u[i] = s.get1d();
    }
	TMarkovChain large_step(sampler &s) const
    {
		TMarkovChain Result;
		Result.C = (*this).C;
		for (int i = 0; i < NumStates; i++) Result.u[i] = s.get1d();
		return Result;
	}
	TMarkovChain mutate(sampler& s) const {
		TMarkovChain Result;
		Result.C = (*this).C;

		// pixel location
		Result.u[0] = perturb(u[0], 2.0 / float(PixelWidth+PixelHeight), 0.1f, s);
		Result.u[1] = perturb(u[1], 2.0 / float(PixelWidth+PixelHeight), 0.1f, s);

		// the rest
		for (int i = 2; i < NumStates; i++) Result.u[i] = perturb(u[i], 1.0f / 1024.0f, 1.0f / 64.0f, s);
		return Result;
	}
};

void InitRandomNumbersByChain(const TMarkovChain &MC, float *prnds)
{ 
    for (int i = 0; i < NumStates; i++) 
        prnds[i] = MC.u[i]; 
}

void InitRandomNumbers(sampler &s, float *prnds) 
{ 
    for (int i = 0; i < NumStates; i++)
        prnds[i] = s.get1d();
}

// path sampling
void pssmlt::TracePath(pssmlt::Path &path, const ray &r, Scene *scene, float *prnds, int &PathRndsOffset)
{
    hit_record dummy_hrec;
    // Compute scalar contribution function for this path
    state st;
    st.depth = 0;
    st.prev_bsdf_pdf = 0.0f;
    st.prnds = prnds;
    st.PathRndsOffset = PathRndsOffset;
    st.prev_hrec = dummy_hrec;
    Vector3f c = Li(path, r, scene, st);
    PathRndsOffset = st.PathRndsOffset;
    path.contrib.c = c;
    path.contrib.sc = std::max(std::max(c[0], c[1]), c[2]);
}

pssmlt::Path pssmlt::GenerateEyePath(const int MaxEyeEvents, Scene *scene, float *prnds, int &PathRndsOffset)
{
    Path Result;
    Result.n = 0;
    camera &cam = scene->cam;

    if (MaxEyeEvents == 0) return Result;
    for (int i = 0; i < MaxEvents; i++) Result.x[i].hit_obj = (hitable*)0xdeadbeef;
    PathRndsOffset = 0;

    ray r = cam.get_ray(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1], sampler().get2d());
    
    //const Vector3f su = cam.u * -(0.5 - prnds[PathRndsOffset + 0]) * PixelWidth;
    //const Vector3f sv = cam.v * (0.5 - prnds[PathRndsOffset + 1]) * PixelHeight;
    //const Vector3f sw = cam.w * cam.dist;
    //ray r = ray(cam.origin, unit_vector(su + sv + sw));
    //std::cout << r.d[0] << " " << r.d[1] << " " << r.d[2] << std::endl;


    PathRndsOffset += 2;

    Result.x[0] = Vert(r.o, cam.w, nullptr); 
    Result.camera_ray = unit_vector(r.d);

    Result.n++;
    hit_record dummy_hrec;
    TracePath(Result, r, scene, prnds, PathRndsOffset);

    // get the pixel location
    Vector3f Direction = Result.camera_ray;
    Vector3f ScreenCenter = cam.origin + (cam.w * cam.dist);
    Vector3f ScreenPosition = cam.origin + (Direction * (cam.dist / dot(Direction, cam.w))) - ScreenCenter;
    Result.contrib.x = -dot(cam.u, ScreenPosition) + (PixelWidth * 0.5f);
    Result.contrib.y = -dot(cam.v, ScreenPosition) + (PixelHeight * 0.5f);

    //std::cout << Direction[0] << " " << Direction[1] << " " << Direction[2] << std::endl;
    //std::cout << Result.contrib.x << " " << Result.contrib.y << std::endl;

    return Result;
}


Vector3f pssmlt::Li(Path &path, const ray &r, Scene *scene, state &st)
{
    hit_record hrec;
    auto &world = scene->world;
    auto &lights = scene->lights;
    if (world->hit(r, EPSILON, FLT_MAX, hrec))
    {
        scatter_record srec;
        Vector3f Le = hrec.mat_ptr->emitted(r, hrec);

        // set path data
        path.x[path.n] = Vert(hrec.p, hrec.normal, hrec.obj); 
        path.n++;
        float rnd0 = st.prnds[st.PathRndsOffset + 0];
        float rnd1 = st.prnds[st.PathRndsOffset + 1];
        float rnd2 = st.prnds[st.PathRndsOffset + 2];
        st.PathRndsOffset += NumRNGsPerEvent;
        Vector3f rnd(rnd0, rnd1, rnd2);

        /* If we hit a light source, weight its contribution */
        if (((Le.r() != 0.0f) || (Le.g() != 0.0f) || (Le.b() != 0.0f)))
        {
            if ((st.depth == 0)
                || (dynamic_cast<modified_phong*>(st.prev_hrec.mat_ptr) != nullptr)
                || (dynamic_cast<metal*>(st.prev_hrec.mat_ptr) != nullptr))
                return Le;
            // Start with checking if camera ray hits a light source
            const float cos_wo = std::max(dot(hrec.normal, -unit_vector(r.direction())), 0.0f);
            float distance_squared = (hrec.p - st.prev_hrec.p).squared_length();
            if (distance_squared <= EPSILON) distance_squared = EPSILON;

            const float surface_bsdf_pdf = st.prev_bsdf_pdf * cos_wo / distance_squared;
            const float light_pdf = lights.pdf_direct_sampling(hrec, r.direction());

            const float weight = miWeight(surface_bsdf_pdf, light_pdf);

            return Le * weight;
        }

        if (st.depth <= MaxPathLength && hrec.mat_ptr->scatter(r, hrec, srec, rnd))
        {
            if (srec.is_specular)
            {
                const float surface_bsdf_pdf = srec.pdf_ptr ? srec.pdf_ptr->value(hrec, srec.specular_ray.direction()) : 1.0f;
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(r, hrec, srec.specular_ray.direction());
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                //const float cos_wi = abs(dot(hrec.normal, unit_vector(srec.specular_ray.direction())));
                srec.specular_ray.o += (EPSILON * hrec.normal);
                st.depth += 1;
                st.prev_bsdf_pdf = surface_bsdf_pdf;
                st.prev_hrec = hrec;
                
                return surface_bsdf * Li(path, srec.specular_ray, scene, st) / surface_bsdf_pdf;
            }
            else
            {
                rnd0 = st.prnds[st.PathRndsOffset + 0];
                rnd1 = st.prnds[st.PathRndsOffset + 1];
                rnd2 = st.prnds[st.PathRndsOffset + 2];
                st.PathRndsOffset += NumRNGsPerEvent;
                /* Direct light sampling */
                const int index = lights.pick_sample(rnd0);
                if (index >= 0)
                {
                    /* Sample a random light source */
                    hit_record lrec;
                    Vector3f offset_origin = hrec.p + (EPSILON * hrec.normal);

                    Vector2f rnd_2d(rnd1, rnd2);
                    Vector3f to_light = lights[index]->sample_direct(lrec, offset_origin, rnd_2d);
                    const float dist_to_light = to_light.length();
                    //to_light.make_unit_vector();

                    ray shadow_ray = ray(offset_origin, to_light);

                    if (!world->hit(shadow_ray, EPSILON, 1 - SHADOW_EPSILON, lrec))
                    {
                        const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(shadow_ray, hrec, to_light);
                        // Calculate geometry term
                        const float cos_wi = std::abs(dot(hrec.normal, unit_vector(to_light)));
                        const float cos_wo = std::max(dot(lrec.normal, -unit_vector(to_light)), 0.0f);
                        if (cos_wo != 0)
                        {
                            float distance_squared = dist_to_light * dist_to_light;

                            const float light_pdf = lights.pdf_direct_sampling(hrec, to_light);
                            // Visibility term is always 1
                            // because of the invariant imposed on these objects by the if above.
                            const float surface_bsdf_pdf = srec.pdf_ptr->value(hrec, to_light) * cos_wo / distance_squared;

                            if (distance_squared <= EPSILON) distance_squared = EPSILON;

                            const float G = cos_wi * cos_wo / distance_squared;

                            const float weight = miWeight(light_pdf, surface_bsdf_pdf);

                            Le += lights.list_size * lrec.mat_ptr->emitted(shadow_ray, lrec) * surface_bsdf * G * weight / light_pdf;
                        }
                    }
                }
                /* Sample BSDF to generate next ray direction for indirect lighting */
                hrec.p = hrec.p + (EPSILON * hrec.normal);

                Vector2f rnd_2d(st.prnds[st.PathRndsOffset + 0], st.prnds[st.PathRndsOffset + 1]);
                st.PathRndsOffset += 2;
                ray wo(hrec.p, srec.pdf_ptr->generate(rnd_2d, hrec));
                const float surface_bsdf_pdf = srec.pdf_ptr->value(hrec, wo.direction());
                const Vector3f surface_bsdf = hrec.mat_ptr->eval_bsdf(wo, hrec, wo.direction());
                /* Reject current path in case the ray is on the wrong side of the surface (BRDF is 0 as ray is pointing away from the hemisphere )*/
                if (surface_bsdf_pdf == 0)
                {
                    return Vector3f(0, 0, 0);
                }
                const float cos_wo = abs(dot(hrec.normal, unit_vector(wo.direction())));
                st.depth += 1;
                st.prev_bsdf_pdf = surface_bsdf_pdf;
                st.prev_hrec = hrec;

                return Le + surface_bsdf * Li(path, wo, scene, st) * cos_wo / surface_bsdf_pdf;
            }
        }
        return Le;
    }

    //Vector3f unit_direction = unit_vector(r.direction());
    //float t = 0.5*(unit_direction.y() + 1.0);
    //return (1.0 - t)*Vector3f(1.0, 1.0, 1.0) + t * Vector3f(0.5, 0.7, 1.0);
    return scene->env_map->eval(r, st.prev_hrec, st.depth);
}

void pssmlt::background_thread(const std::shared_future<void> &future, GLFWwindow *window, viewer *film_viewer)
{
    while (!film_viewer->to_exit)
    {
        /* Poll for and process events */
        glfwPollEvents();
        /* Render here */
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, film_viewer->nx, film_viewer->ny, 0, GL_RGB, GL_UNSIGNED_BYTE, film_viewer->out_image.get());
        glBlitFramebuffer(0, 0, film_viewer->nx, film_viewer->ny, 0, 0, film_viewer->nx, film_viewer->ny,
            GL_COLOR_BUFFER_BIT, GL_NEAREST);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        if (future.wait_for(std::chrono::seconds(0)) == std::future_status::ready || glfwWindowShouldClose(window))
        {
            film_viewer->to_exit = true;
        }
    }
}


void pssmlt::Render(Scene *scene, viewer *film_viewer, tf::Taskflow &tf)
{
    sampler init_s;
    double b = 0.0;
    float prnds_init[NumStates];
    int PathRndsOffsetInit;
    for (int i = 0; i < N_Init; i++)
    {
        InitRandomNumbers(init_s, prnds_init);
        b += GenerateEyePath(MaxEvents, scene, prnds_init, PathRndsOffsetInit).contrib.sc;
    }
    b /= N_Init;

    //fprintf(stderr, "\n");

    tf.parallel_for((int)tf.num_workers(), 0, -1, [&](int thread_id)
    {
        // internal states of random numbers
        int PathRndsOffset;
        float prnds[NumStates];

        static thread_local sampler s(thread_id * 39);
        std::cout << thread_id << std::endl;
        TMarkovChain current(s), proposal(s);
        InitRandomNumbersByChain(current, prnds);
        current.C = GenerateEyePath(MaxEvents, scene, prnds, PathRndsOffset).contrib;
        //fprintf(stderr, "\r%f%%", (float(total_samples) / film_viewer->ns)*100.0f);
        const int samples_per_thread = film_viewer->ns / tf.num_workers();
        for (int total_samples = 0; total_samples < samples_per_thread; ++total_samples)
        {
            if (film_viewer->to_exit) break;
            float is_large_step_done;
            if (s.get1d() < LargeStepProb)
            {
                proposal = current.large_step(s);
                is_large_step_done = 1.0f;
            }
            else
            {
                proposal = current.mutate(s);
                is_large_step_done = 0.0f;
            }
            InitRandomNumbersByChain(proposal, prnds);
            proposal.C = GenerateEyePath(MaxEvents, scene, prnds, PathRndsOffset).contrib;

            float a = 1.0f;
            if (current.C.sc > 0.0)
                a = std::max(std::min(1.0, proposal.C.sc / current.C.sc), 0.0);

            // accumulate samples
            if (proposal.C.sc > 0.0) AccumulatePathContribution(proposal.C, (a + is_large_step_done) / (proposal.C.sc / b + LargeStepProb), film_viewer);
            if (current.C.sc > 0.0) AccumulatePathContribution(current.C, (1.0 - a) / (current.C.sc / b + LargeStepProb), film_viewer);

            // update the chain
            if (s.get1d() <= a) current = proposal;

            //std::shared_future<void> dummy_future;
            //if (total_samples % 100000 == 0)
            //    background_thread(dummy_future, film_viewer->window, film_viewer);
        }
    });
    auto future = tf.dispatch();
    background_thread(future, film_viewer->window, film_viewer);
    std::cout << "\ndone";
}