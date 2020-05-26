#pragma once

#include <memory>
#include <vector>
#include <cassert>
#include <cstring>
#include "hitable.h"
#include "material.h"
#include "SHSample.h"
#include "Scene.h"
#include "triangle.h"

class triangle_sse : public triangle
{
public:
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const { std::cout << "hit Not implemented in triangle_sse!"; return false; }
    //Use Moller-Trumbore intersection algorithm (Fast Minimum Storage Ray/Triangle Intersection)
    virtual int hit_simd(ray4 &r4, float t_min, float t_max, hit_record4 &hrec4) const override
    {
        ++Scene::num_ray_tri_intersections;
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];

        // Vector3f e1 = v1 - v0;
        // Vector3f e2 = v2 - v0;
        __m128 e1x4 = _mm_set1_ps(v1.x() - v0.x());
        __m128 e1y4 = _mm_set1_ps(v1.y() - v0.y());
        __m128 e1z4 = _mm_set1_ps(v1.z() - v0.z());
        __m128 e2x4 = _mm_set1_ps(v2.x() - v0.x());
        __m128 e2y4 = _mm_set1_ps(v2.y() - v0.y());
        __m128 e2z4 = _mm_set1_ps(v2.z() - v0.z());

        // const Vector3f h = cross(r.d, edge2);
        __m128 Px4 = _mm_sub_ps(
            _mm_mul_ps(r4.dy4, e2z4),
            _mm_mul_ps(r4.dz4, e2y4)
        );
        __m128 Py4 = _mm_sub_ps(
            _mm_mul_ps(r4.dz4, e2x4),
            _mm_mul_ps(r4.dx4, e2z4)
        );
        __m128 Pz4 = _mm_sub_ps(
            _mm_mul_ps(r4.dx4, e2y4),
            _mm_mul_ps(r4.dy4, e2x4)
        );

        //a = dot(edge1, h);
        __m128 det4 = _mm_add_ps(
            _mm_add_ps(
                _mm_mul_ps(e1x4, Px4),
                _mm_mul_ps(e1y4, Py4)
            ),
            _mm_mul_ps(e1z4, Pz4)
        );

        // Check if this ray is parallel to this triangle's plane.
        //if (a == 0)
        //    return false;
        __m128 mask1 = _mm_or_ps(
            _mm_cmple_ps(det4, MINUSEPS4),
            _mm_cmpge_ps(det4, EPS4)
        );

        // f = 1.0f / a;
        __m128 inv_det4 = _mm_rcp_ps(det4);

        // Vector3f s = r.o - v0;
        __m128 Tx4 = _mm_sub_ps(r4.ox4, _mm_set_ps1(v1.x));
        __m128 Ty4 = _mm_sub_ps(r4.oy4, _mm_set_ps1(v1.y));
        __m128 Tz4 = _mm_sub_ps(r4.oz4, _mm_set_ps1(v1.z));
        __m128 u4 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(Tx4, Px4), _mm_mul_ps(Ty4, Py4)), _mm_mul_ps(Tz4, Pz4)), inv_det4);
        __m128 mask2 = _mm_and_ps(_mm_cmpge_ps(u4, _mm_setzero_ps()), _mm_cmple_ps(u4, ONE4));
        __m128 Qx4 = _mm_sub_ps(_mm_mul_ps(Ty4, e1z4), _mm_mul_ps(Tz4, e1y4));
        __m128 Qy4 = _mm_sub_ps(_mm_mul_ps(Tz4, e1x4), _mm_mul_ps(Tx4, e1z4));
        __m128 Qz4 = _mm_sub_ps(_mm_mul_ps(Tx4, e1y4), _mm_mul_ps(Ty4, e1x4));
        __m128 v4 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(r4.dx4, Qx4), _mm_mul_ps(r4.dy4, Qy4)), _mm_mul_ps(r4.dz4, Qz4)), inv_det4);
        __m128 mask3 = _mm_and_ps(_mm_cmpge_ps(v4, _mm_setzero_ps()), _mm_cmple_ps(_mm_add_ps(u4, v4), ONE4));
        __m128 t4 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(e2x4, Qx4), _mm_mul_ps(e2y4, Qy4)), _mm_mul_ps(e2z4, Qz4)), inv_det4);
        __m128 mask4 = _mm_cmpgt_ps(t4, _mm_setzero_ps());
        __m128 mask5 = _mm_cmplt_ps(t4, r4.t4);

        //__m128 mask6 = _mm_and_ps(
        //    _mm_cmpge_ps(det4, _mm_set_ps1(t_min)),
        //    _mm_cmple_ps(det4, _mm_set_ps1(t_max))
        //);

        __m128 combined_mask = _mm_and_ps(_mm_and_ps(_mm_and_ps(_mm_and_ps(mask1, mask2), mask3), mask4), mask5);

        // hrec.t = t;
        r4.t4 = _mm_blendv_ps(r4.t4, t4, combined_mask);
        hrec4.t4 = r4.t4;

        // hrec.wi = -unit_vector(r.d);
        __m128 sqX4 = _mm_mul_ps(r4.ox4, r4.ox4);
        __m128 sqY4 = _mm_mul_ps(r4.oy4, r4.oy4);
        __m128 sqZ4 = _mm_mul_ps(r4.oz4, r4.oz4);
        __m128 len4 = _mm_rsqrt_ps(_mm_add_ps(_mm_add_ps(sqX4, sqY4), sqZ4)); // reciprocal square root

        // A lot of unreadable code, but it basically combines negation and multiplication with reciprocal sqrt
        hrec4.wix4 = _mm_xor_ps(_mm_mul_ps(r4.ox4, len4), _mm_set1_ps(-0.f));
        hrec4.wiy4 = _mm_xor_ps(_mm_mul_ps(r4.oy4, len4), _mm_set1_ps(-0.f));
        hrec4.wiz4 = _mm_xor_ps(_mm_mul_ps(r4.oz4, len4), _mm_set1_ps(-0.f));

        // Use u, v to find interpolated normals and texture coords
        // P = (1 - u - v) * V0 + u * V1 + v * V2
        //hrec.normal = unit_vector(cross(edge1, edge2));
        //hrec.normal = unit_vector((1 - u - v) * mesh->normals[V[0]] + u * mesh->normals[V[1]] + v * mesh->normals[V[2]]);
        const Vector3f &n0 = mesh->normals[V[0]];
        const Vector3f &n1 = mesh->normals[V[1]];
        const Vector3f &n2 = mesh->normals[V[2]];

        __m128 n0x4 = _mm_set1_ps(n0.x());
        __m128 n0y4 = _mm_set1_ps(n0.y());
        __m128 n0z4 = _mm_set1_ps(n0.z());
        __m128 n1x4 = _mm_set1_ps(n1.x());
        __m128 n1y4 = _mm_set1_ps(n1.y());
        __m128 n1z4 = _mm_set1_ps(n1.z());
        __m128 n2x4 = _mm_set1_ps(n2.x());
        __m128 n2y4 = _mm_set1_ps(n2.y());
        __m128 n2z4 = _mm_set1_ps(n2.z());
        
        // (1-u-v) factor i.e., w for the u,v
        __m128 w4 = _mm_sub_ps(_mm_set_ps1(1.0f), _mm_sub_ps(u4, v4));

        hrec4.nx4 = _mm_add_ps(_mm_mul_ps(w4, n0x4), _mm_add_ps(_mm_mul_ps(u4, n1x4), _mm_mul_ps(v4, n2x4)));
        hrec4.ny4 = _mm_add_ps(_mm_mul_ps(w4, n0y4), _mm_add_ps(_mm_mul_ps(u4, n1y4), _mm_mul_ps(v4, n2y4)));
        hrec4.nz4 = _mm_add_ps(_mm_mul_ps(w4, n0z4), _mm_add_ps(_mm_mul_ps(u4, n1z4), _mm_mul_ps(v4, n2z4)));

        //Vector2f uvhit = (1 - u - v) * mesh->uv[V[0]] + u * mesh->uv[V[1]] + v * mesh->uv[V[2]];
        //hrec.u = uvhit.x;
        //hrec.v = uvhit.y;
        const Vector2f &uv0 = mesh->uv[V[0]];
        const Vector2f &uv1 = mesh->uv[V[1]];
        const Vector2f &uv2 = mesh->uv[V[2]];

        __m128 u0_4 = _mm_set1_ps(uv0.x);
        __m128 v0_4 = _mm_set1_ps(uv0.y);
        __m128 u1_4 = _mm_set1_ps(uv1.x);
        __m128 v1_4 = _mm_set1_ps(uv1.y);
        __m128 u2_4 = _mm_set1_ps(uv2.x);
        __m128 v2_4 = _mm_set1_ps(uv2.y);

        hrec4.u4 = _mm_add_ps(_mm_mul_ps(w4, u0_4), _mm_add_ps(_mm_mul_ps(u4, u1_4), _mm_mul_ps(v4, u2_4)));
        hrec4.v4 = _mm_add_ps(_mm_mul_ps(w4, v0_4), _mm_add_ps(_mm_mul_ps(u4, v1_4), _mm_mul_ps(v4, v2_4)));

        hrec4.mat_ptr = mat_ptr;
        hrec4.obj = (hitable *)this;


        return _mm_movemask_ps(combined_mask);
    }

    virtual bool bounding_box(float t0, float t1, aabb &b) const
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

    virtual float pdf_direct_sampling(const hit_record &lrec, const Vector3f &to_light) const
    {
        // This is explicitly converting to area measure
        // TODO: Allow switching between solid angle/area measure
        return inv_area;
    }
    virtual Vector3f sample_direct(hit_record &rec, const Vector3f &o, const Vector2f &sample) const
    {
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];

        const Vector2f &u = sample;
        float su0 = std::sqrt(u.x);
        float b0 = 1 - su0;
        float b1 = u.y * su0;

        Vector3f random_point = (1 - b0 - b1) * v0 + b0 * v1 + b1 * v2;

        rec.t = 1.0f;
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

    // Triangle vertex index data
    const int *V;
    float inv_area;
    // Store edges here so we don't calculate for every intersection test
    const Vector3f edge1;
    const Vector3f edge2;
    std::shared_ptr<triangle_mesh> mesh;
    material *mat_ptr;
    // For PRT
    std::unique_ptr<PRT::SHCoefficients[]> coeffs;
};