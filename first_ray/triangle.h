#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <memory>
#include <vector>
#include <cassert>
#include <cstring>
#include "hitable.h"
#include "material.h"
#include "SHSample.h"
#include "Scene.h"


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
    }

    //Use M�ller�Trumbore intersection algorithm (Fast Minimum Storage Ray/Triangle Intersection)
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &hrec) const
    {
        ++Scene::num_ray_tri_intersections;
        double a, f, u, v;
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f h = cross(r.d, edge2);
        a = dot(edge1, h);
        // Check if this ray is parallel to this triangle's plane.
        if (a == 0)
            return false;    
        f = 1.0 / a;
        Vector3f s = r.o - v0;
        u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        Vector3f q = cross(s, edge1);
        v = f * dot(r.d, q);
        if (v >= 0.0 && u + v <= 1.0)
        {
            // At this stage we can compute t to find out where the intersection point is on the line.
            double t = f * dot(edge2, q);

            // Check if ray intersected successfully
            if (t > t_min
                && t < t_max)
            {
                hrec.t = t;
                hrec.p = r.point_at_parameter(t);
                // Use u, v to find interpolated normals and texture coords
                // P = (1 - u - v) * V0 + u * V1 + v * V2
                //hrec.normal = unit_vector(cross(edge1, edge2));
                if (mesh->use_geometry_normals)
                    hrec.normal = unit_vector(cross(edge1, edge2));
                else
				    hrec.normal = unit_vector((1 - u - v) * mesh->normals[V[0]] + u * mesh->normals[V[1]] + v * mesh->normals[V[2]]);

                Vector2f uvhit = (1 - u - v) * mesh->uv[V[0]] + u * mesh->uv[V[1]] + v * mesh->uv[V[2]];
                hrec.u = uvhit.x;
                hrec.v = uvhit.y;
                hrec.wi = -unit_vector(r.d);
				hrec.uv.x = u;
				hrec.uv.y = v;
                hrec.mat_ptr = mat_ptr;
				hrec.obj = (hitable*)this;
                return true;
            }
        }
        // This means that there is a line intersection but not a ray intersection.
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
        if (mesh->use_geometry_normals)
            rec.normal = unit_vector(cross(edge1, edge2));
        else
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
std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, const Matrix4x4 &toWorld, material *bsdf, std::vector<hitable *> &lights, bool use_geometry_normals = false);

#endif