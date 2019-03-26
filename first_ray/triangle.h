#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "hitable.h"
#include <memory>
#include <vector>

class triangle_mesh
{
public:

    triangle_mesh(const int &nTriangles, const int &nVertices, Vector3f *v,
        const int *indices, Vector3f *n, Point2f *_uv, const bool &shallow_copy = false)
        : nTriangles(nTriangles),
        nVertices(nVertices),
        indices(indices, indices + 3 * nTriangles)

    {

        vertices.reset(shallow_copy ? v : new Vector3f[nVertices]);
        normals.reset(shallow_copy ? n : new Vector3f[nVertices]);
        uv.reset(shallow_copy ? _uv : new Point2f[nVertices]);

        if (!shallow_copy)
        {
            memcpy(vertices.get(), v, sizeof(Vector3f) * nVertices);
            memcpy(normals.get(), v, sizeof(Vector3f) * nVertices);
            memcpy(uv.get(), v, sizeof(Point2f) * nVertices);
        }

    }

    int nTriangles;
    int nVertices;
    std::vector<int> indices;
    std::unique_ptr<Vector3f[]> vertices;
    std::unique_ptr<Vector3f[]> normals;
    std::unique_ptr<Point2f[]> uv;
};

class triangle : public hitable
{
public:
    triangle(const std::shared_ptr<triangle_mesh> &mesh, const int &tri_num, std::shared_ptr<material> &mat) : mesh(mesh), mat_ptr(mat)
    {
        V = &mesh->indices[3 * tri_num];
    }
    //triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2, material *mat) : v0(_v0), v1(_v1), v2(_v2), mat_ptr(mat) {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        const float EPSILON = 0.0000001;
        float a, f, u, v;
        const Vector3f v0 = mesh->vertices[V[0]];
        const Vector3f v1 = mesh->vertices[V[1]];
        const Vector3f v2 = mesh->vertices[V[2]];
        const Vector3f edge1 = v1 - v0;
        const Vector3f edge2 = v2 - v0;
        const Vector3f h = cross(r.d, edge2);
        a = dot(edge1, h);
        // Check if this ray is parallel to this triangle's plane.
        if (a > -EPSILON && a < EPSILON)
            return false;    
        f = 1.0 / a;
        Vector3f s = r.o - v0;
        u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        Vector3f q = cross(s, edge1);
        v = f * dot(r.d, q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        float t = f * dot(edge2, q);

        // Check if ray intersected successfully
        if (t > EPSILON)
        {
            rec.t = t;
            rec.p = r.point_at_parameter(t);
            rec.normal = cross(edge1, edge2);
            rec.normal.make_unit_vector();
            rec.u = u;
            rec.v = v;
            rec.mat_ptr = mat_ptr.get();
            return true;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return false;
    }

    // Triangle data
    const int *V;
    std::shared_ptr<triangle_mesh> mesh;

    std::shared_ptr<material> mat_ptr;
};

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, std::shared_ptr<material> mat);


#endif