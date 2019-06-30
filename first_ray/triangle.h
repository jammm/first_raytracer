#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <memory>
#include <vector>
#include <cassert>
#include <cstring>
#include "hitable.h"
#include "material.h"


class triangle_mesh
{
public:

	triangle_mesh() {}

    triangle_mesh(const int &nTriangles, const int &nVertices, Vector3f *v,
        const int *indices, Vector3f *n, Point2f *_uv, std::unique_ptr<material> mat, std::string name, const bool &shallow_copy = false)
        : nTriangles(nTriangles),
        nVertices(nVertices),
        indices(indices, indices + 3 * nTriangles), mat(move(mat)),
        name(name)
    {

        vertices.reset(shallow_copy ? v : new Vector3f[nVertices]);
        normals.reset(shallow_copy ? n : new Vector3f[nVertices]);
        uv.reset(shallow_copy ? _uv : new Point2f[nVertices]);

        if (!shallow_copy)
        {
            std::memcpy(vertices.get(), v, sizeof(Vector3f) * nVertices);
            std::memcpy(normals.get(), v, sizeof(Vector3f) * nVertices);
            std::memcpy(uv.get(), v, sizeof(Point2f) * nVertices);
        }

    }

    int nTriangles;
    int nVertices;
    std::vector<int> indices;
    std::unique_ptr<Vector3f[]> vertices;
    std::unique_ptr<Vector3f[]> normals;
    std::unique_ptr<Point2f[]> uv;
    std::unique_ptr<material> mat;
    // Name of the object this mesh belongs to
    std::string name;
};

class triangle : public hitable
{
public:
    triangle(const std::shared_ptr<triangle_mesh> &mesh, const int &tri_num) : mesh(mesh),
        mat_ptr(mesh->mat.get()),
        edge1(mesh->vertices[mesh->indices[3 * tri_num + 1]] - mesh->vertices[mesh->indices[3 * tri_num]]),
        edge2(mesh->vertices[mesh->indices[3 * tri_num + 2]] - mesh->vertices[mesh->indices[3 * tri_num]])
    {
        V = &mesh->indices[3 * tri_num];
    }

    //Use M�ller�Trumbore intersection algorithm (Fast Minimum Storage Ray/Triangle Intersection)
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        float a, f, u, v;
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];
        const Vector3f h = cross(r.d, edge2);
        a = dot(edge1, h);
        // Check if this ray is parallel to this triangle's plane.
        if (a > -EPSILON && a < EPSILON)
            return false;    
        f = 1.0f / a;
        Vector3f s = r.o - v0;
        u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        Vector3f q = cross(s, edge1);
        v = f * dot(r.d, q);
        if (v >= 0.0f && u + v <= 1.0)
        {
            // At this stage we can compute t to find out where the intersection point is on the line.
            float t = f * dot(edge2, q);

            // Check if ray intersected successfully
            if (t > t_min
                && t < t_max)
            {
                rec.t = t;
                rec.p = r.point_at_parameter(t);
                rec.normal = cross(edge1, edge2);
                rec.normal.make_unit_vector();
                // Use u, v to find interpolated coordinates
                // P = (1 - u - v) * V0 + u * V1 + v * V2 
                Point2f uvhit = (1 - u - v) * mesh->uv[V[0]] + u * mesh->uv[V[1]] + v * mesh->uv[V[2]];
                rec.u = uvhit.x;
                rec.v = uvhit.y;
                rec.mat_ptr = mat_ptr.get();
                return true;
            }
        }
        // This means that there is a line intersection but not a ray intersection.
        return false;
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
        const float area = 0.5f * cross(edge1, edge2).length();

        // This is explicitly converting to area measure
        // TODO: Allow switching between solid angle/area measure
        return 1.0f / area;
    }
    virtual Vector3f random(const Vector3f &o) const
    {
        const Vector3f &v0 = mesh->vertices[V[0]];
        const Vector3f &v1 = mesh->vertices[V[1]];
        const Vector3f &v2 = mesh->vertices[V[2]];

        Point2f u(gen_cano_rand(), gen_cano_rand());
        float su0 = std::sqrt(u.x);
        float b0 = 1 - su0;
        float b1 = u.y * su0;

        Vector3f random_point = v0 * b0 + v1 * b1 + v2 * (1 - b0 - b1);

        return random_point - o;
    }

    // Triangle vertex index data
    const int *V;
    // Store edges here so we don't calculate for every intersection test
    const Vector3f edge1;
    const Vector3f edge2;
    std::shared_ptr<triangle_mesh> mesh;
    std::shared_ptr<material> mat_ptr;
};

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, std::vector<hitable *> &lights);


#endif