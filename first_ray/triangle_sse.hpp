#pragma once

#include <memory>
#include <vector>
#include <cassert>
#include <cstring>
#include "hitable.h"
#include "material.h"
#include "SHSample.h"
#include "Scene.h"

#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__ ((aligned (16)))
#endif

// Yet Faster Ray-Triangle Intersection (Using SSE4)
// Jiri Havel and Adam Herout

class triangle_mesh
{
public:

	triangle_mesh() {}

    triangle_mesh(const int &nTriangles, const int &nVertices, Vector3f *v,
        const int *indices, Vector3f *n, Vector2f *_uv, std::unique_ptr<material> mat, bool _use_geometry_normals,
        std::string name, const bool &shallow_copy = false)
        : nTriangles(nTriangles),
        nVertices(nVertices),
        indices(indices, indices + 3 * nTriangles), mat(move(mat)),
        use_geometry_normals(_use_geometry_normals),
        name(name)
    {

        vertices.reset(shallow_copy ? v : new Vector3f[nVertices]);
        normals.reset(shallow_copy ? n : new Vector3f[nVertices]);
        uv.reset(shallow_copy ? _uv : new Vector2f[nVertices]);

        if (!shallow_copy)
        {
            std::memcpy(vertices.get(), v, sizeof(Vector3f) * nVertices);
            std::memcpy(normals.get(), v, sizeof(Vector3f) * nVertices);
            std::memcpy(uv.get(), v, sizeof(Vector2f) * nVertices);
        }

    }

    int nTriangles;
    int nVertices;
    std::vector<int> indices;
    std::unique_ptr<Vector3f[]> vertices;
    std::unique_ptr<Vector3f[]> normals;
    std::unique_ptr<Vector2f[]> uv;
    std::unique_ptr<material> mat;
    bool use_geometry_normals;
    // Name of the object this mesh belongs to
    std::string name;
};

struct alignas(16) PrecomputedTriangle
{
    double nx, ny, nz, nd;
    double ux, uy, uz, ud;
    double vx, vy, vz, vd;
};

const double int_coef_arr[4] = { -1, -1, -1, 1 };
const __m128 int_coef = _mm_load_ps(int_coef_arr);

const double inv_coef_p_arr[4] = { 1, 1, 1, 1 };
const __m128 inv_coef_p = _mm_load_ps(inv_coef_p_arr);

class triangle : public hitable
{
public:
    triangle(const std::shared_ptr<triangle_mesh> &mesh, const int &tri_num) :
        edge1(mesh->vertices[mesh->indices[3 * tri_num + 1]] - mesh->vertices[mesh->indices[3 * tri_num]]),
        edge2(mesh->vertices[mesh->indices[3 * tri_num + 2]] - mesh->vertices[mesh->indices[3 * tri_num]]),
        mesh(mesh),
        mat_ptr(mesh->mat.get())
    {
        V = &mesh->indices[3 * tri_num];
        inv_area = 1 / (0.5 * cross(edge1, edge2).length() * mesh->nTriangles);

        Vector3f A = mesh->vertices[V[0]];

        Vector3f n = cross(edge1, edge2);
        if (n.squared_length() == 0)
            n = mesh->normals[V[0]];

        double n_sq_len = n.squared_length();

        p.nx = n.x();
        p.ny = n.y();
        p.nz = n.z();
        p.nd = dot(A, n);

        const double inv_n_sq_len = 1.0 / n_sq_len;

        Vector3f n1 = cross(edge2, n) * inv_n_sq_len;
        p.ux = n1.x();
        p.uy = n1.y();
        p.uz = n1.z();
        p.ud = -dot(n1, A);

        Vector3f n2 = cross(n, edge1) * inv_n_sq_len;
        p.vx = n2.x();
        p.vy = n2.y();
        p.vz = n2.z();
        p.vd = -dot(n2, A);
    }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &hrec) const
    {
        ++Scene::num_ray_tri_intersections;
        ray r_new = r;
        r_new.o.e[3] = 1.0;
        r_new.d.e[3] = 0.0;
        const __m128 o = _mm_load_ps(&(r_new.o.e[0]));
        const __m128 d = _mm_load_ps(&(r_new.d.e[0]));
        const __m128 n = _mm_load_ps(&(p.nx));
        const __m128 det = _mm_dp_ps(n, d, 0x7f);
        const __m128 dett = _mm_dp_ps(
            _mm_mul_ps(int_coef, n), o, 0xff);
        {
            const __m128 detp = _mm_add_ps(_mm_mul_ps(o, det),
                _mm_mul_ps(dett, d));
            const __m128 detu = _mm_dp_ps(detp,
                _mm_load_ps(&p.ux), 0xf1);
            if ((_mm_movemask_ps(_mm_xor_ps(detu,
                _mm_sub_ss(det, detu))) & 1) == 0)
            {
                const __m128 detv = _mm_dp_ps(detp, _mm_load_ps(&p.vx), 0xf1);
                if ((_mm_movemask_ps(_mm_xor_ps(detv,
                    _mm_sub_ss(det, _mm_add_ss(detu, detv)))) & 1) == 0)
                {
                    const __m128 inv_det = _mm_rcp_ps(det);
                    _mm_store_ss(&hrec.t, _mm_mul_ss(dett, inv_det));
                    if ((hrec.t > t_min) && (hrec.t < t_max))
                    {
                        _mm_store_ss(&hrec.uv.x, _mm_mul_ss(detu, inv_det));
                        _mm_store_ss(&hrec.uv.y, _mm_mul_ss(detv, inv_det));
                        _mm_store_ps(&hrec.p.e[0], _mm_mul_ps(detp,
                            _mm_shuffle_ps(inv_det, inv_det, 0)));

                        const double u = hrec.uv.x;
                        const double v = hrec.uv.y;
                        hrec.normal = unit_vector((1 - u - v) * mesh->normals[V[0]] + u * mesh->normals[V[1]] + v * mesh->normals[V[2]]);
                        Vector2f uvhit = (1 - u - v) * mesh->uv[V[0]] + u * mesh->uv[V[1]] + v * mesh->uv[V[2]];
                        hrec.u = uvhit.x;
                        hrec.v = uvhit.y;
                        hrec.wi = -unit_vector(r_new.d);
                        hrec.mat_ptr = mat_ptr;
                        hrec.obj = (hitable *)this;

                        return true;
                    }
                }
            }
        }
        return false;
    }

    virtual bool bounding_box(double t0, double t1, aabb &b) const
    {
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];

        Vector3f minv(fmin(fmin(v0.x(), v1.x()), v2.x()),
                      fmin(fmin(v0.y(), v1.y()), v2.y()),
                      fmin(fmin(v0.z(), v1.z()), v2.z()));

        Vector3f maxv(fmax(fmax(v0.x(), v1.x()), v2.x()),
                      fmax(fmax(v0.y(), v1.y()), v2.y()),
                      fmax(fmax(v0.z(), v1.z()), v2.z()));

        b = aabb(minv, maxv);

        return true;
    }

    virtual double pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const
    {
        // This is explicitly converting to area measure
        // TODO: Allow switching between solid angle/area measure
        return inv_area;
    }
    virtual Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f& sample) const
    {
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];

        const Vector2f &u = sample;
        double su0 = std::sqrt(u.x);
        double b0 = 1 - su0;
        double b1 = u.y * su0;

        Vector3f random_point = (1 - b0 - b1) * v0 + b0 * v1 + b1 * v2;

        rec.t = 1.0;
        rec.p = random_point;
        rec.normal = unit_vector((1 - b0 - b1) * mesh->normals[V[0]] + b0 * mesh->normals[V[1]] + b1 * mesh->normals[V[2]]);
        rec.mat_ptr = mat_ptr;
        Vector2f uvhit = (1 - b0 - b1) * mesh->uv[V[0]] + b0 * mesh->uv[V[1]] + b1 * mesh->uv[V[2]];
        rec.u = uvhit.x;
        rec.v = uvhit.y;
        rec.uv.x = b0;
        rec.uv.y = b1;
		rec.obj = (hitable *)this;
        rec.wi = unit_vector(o - random_point);

        return random_point - o;
    }

    alignas(16) PrecomputedTriangle p;

    // Triangle vertex index data
    const int *V;
    double inv_area;
    // Store edges here so we don't calculate for every intersection test
    const Vector3f edge1;
    const Vector3f edge2;
    std::shared_ptr<triangle_mesh> mesh;
    material *mat_ptr;
	// For PRT
	std::unique_ptr<PRT::SHCoefficients[]> coeffs;
};

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, std::vector<hitable *> &lights, bool use_geometry_normals = false);
std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, const Matrix4x4 &toWorld, material *bsdf, std::vector<hitable*> &lights, bool use_geometry_normals = false);

